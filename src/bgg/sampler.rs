use super::{BggEncoding, BggPublicKey};
use crate::poly::{matrix::*, sampler::*, *};
use itertools::Itertools;
use rand::RngCore;
use std::{marker::PhantomData, sync::Arc};

/// A sampler of a public key A in the BGG+ RLWE encoding scheme
#[derive(Clone)]
pub struct BGGPublicKeySampler<K: AsRef<[u8]>, S: PolyHashSampler<K>> {
    pub params: Arc<<<S::M as PolyMatrix>::P as Poly>::Params>,
    pub sampler: Arc<S>,
    _k: PhantomData<K>,
}

impl<K: AsRef<[u8]>, S: PolyHashSampler<K>> BGGPublicKeySampler<K, S> {
    /// Create a new public key sampler
    /// # Arguments
    /// * `params`: The parameters of the public key matrix
    /// * `sampler`: The sampler to generate the public key matrix
    /// # Returns
    /// A new public key sampler
    pub fn new(params: Arc<<<S::M as PolyMatrix>::P as Poly>::Params>, sampler: Arc<S>) -> Self {
        Self { params, sampler, _k: PhantomData }
    }

    /// Sample a public key matrix
    /// # Arguments
    /// * `tag`: The tag to sample the public key matrix
    /// * `packed_input_size`: The packed input size, i.e., the number of necessary polynomials when
    ///   n bits are encoded into a single polynomial
    /// # Returns
    /// A vector of public key matrices
    pub fn sample(
        &self,
        tag: &[u8],
        packed_input_size: usize,
    ) -> Vec<BggPublicKey<<S as PolyHashSampler<K>>::M>> {
        let log_q = self.params.modulus_bits();
        let columns = 2 * log_q * packed_input_size;
        let all_matrix = self.sampler.sample_hash(tag, 2, columns, DistType::FinRingDist);
        (0..packed_input_size)
            .map(|idx| {
                let matrix = all_matrix.slice_rows(2 * log_q * idx, 2 * log_q * (idx + 1));
                BggPublicKey { matrix }
            })
            .collect_vec()
    }
}

/// A sampler of an encoding in the BGG+ RLWE encoding scheme
///
/// # Fields
/// * `secret`: The secret vector.
/// * `error_sampler`: The sampler to generate LWE errors.
#[derive(Clone)]
pub struct BGGEncodingSampler<S: PolyUniformSampler> {
    pub(crate) secret_vec: S::M,
    pub params: Arc<<<S::M as PolyMatrix>::P as Poly>::Params>,
    pub error_sampler: Arc<S>,
    pub gauss_sigma: f64,
}

impl<S: PolyUniformSampler> BGGEncodingSampler<S> {
    /// Create a new encoding sampler
    /// # Arguments
    /// * `secret`: The secret polynomial
    /// * `params`: The parameters of the encoding
    /// * `error_sampler`: The sampler to generate LWE errors
    /// * `gauss_sigma`: The standard deviation of the Gaussian distribution
    /// # Returns
    /// A new encoding sampler
    pub fn new(
        secret: &<S::M as PolyMatrix>::P,
        params: Arc<<<S::M as PolyMatrix>::P as Poly>::Params>,
        error_sampler: Arc<S>,
        gauss_sigma: f64,
    ) -> Self {
        let minus_one_poly = <S::M as PolyMatrix>::P::const_minus_one(params.as_ref());
        // 1*2 row vector
        let secret_vec = S::M::from_poly_vec_row(&params, vec![secret.clone(), minus_one_poly]);
        Self { secret_vec, params, error_sampler, gauss_sigma }
    }

    pub fn sample(
        &self,
        public_keys: &[BggPublicKey<S::M>],
        plaintexts: &[<S::M as PolyMatrix>::P],
        reveal_plaintexts: bool,
    ) -> Vec<BggEncoding<S::M>> {
        let log_q = self.params.modulus_bits();
        let packed_input_size = plaintexts.len();
        let columns = 2 * log_q * packed_input_size;
        let error = self.error_sampler.sample_uniform(
            1,
            columns,
            DistType::GaussDist { sigma: self.gauss_sigma },
        );
        // first term sA
        // [TODO] Avoid memory cloning here.
        let all_public_key_matrix: S::M = public_keys[0]
            .matrix
            .concat_columns(&public_keys[1..].iter().map(|pk| pk.matrix.clone()).collect_vec());
        let first_term = self.secret_vec.clone() * all_public_key_matrix;
        // second term x \tensor sG
        let gadget = S::M::gadget_matrix(&self.params, 2);
        let sg = self.secret_vec.clone() * gadget;
        let encoded_polys_vec = S::M::from_poly_vec_row(self.params.as_ref(), plaintexts.to_vec());
        let second_term = encoded_polys_vec.tensor(&sg);

        // all_vector = sA + x \tensor sG + e
        let all_vector = first_term + second_term + error;

        let encoding = plaintexts
            .iter()
            .enumerate()
            .map(|(idx, plaintext)| {
                let vector = all_vector.slice_columns(2 * log_q * idx, 2 * log_q * (idx + 1));
                BggEncoding {
                    vector,
                    pubkey: public_keys[idx].clone(),
                    plaintext: if reveal_plaintexts { Some(plaintext.clone()) } else { None },
                }
            })
            .collect_vec();
        encoding
    }
}
