use super::{BggEncoding, BggError, BggPublicKey};
use crate::poly::{gadget::PolyGadgetOps, matrix::*, sampler::*, *};
use itertools::Itertools;
use rand::RngCore;
use std::{marker::PhantomData, sync::Arc};

/// A sampler of a public key A in the BGG+ RLWE encoding scheme
#[derive(Debug, Clone)]
pub struct BGGPublicKeySampler<
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    S: PolyHashSampler<T, P, M, FinRingDist>,
> {
    pub sampler: Arc<S>,
    pub matrix_op: Arc<M>,
    _t: PhantomData<T>,
    _p: PhantomData<P>,
}

impl<
        T: PolyElemOps,
        P: PolyOps<T>,
        M: PolyMatrixOps<T, P>,
        S: PolyHashSampler<T, P, M, FinRingDist>,
    > BGGPublicKeySampler<T, P, M, S>
{
    /// Create a new public key sampler
    /// # Arguments
    /// * `sampler`: The sampler to generate the public key matrix
    /// * `matrix_op`: The matrix operations
    /// # Returns
    /// A new public key sampler
    pub fn new(sampler: Arc<S>, matrix_op: Arc<M>) -> Self {
        Self { sampler, matrix_op, _t: PhantomData, _p: PhantomData }
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
    ) -> Result<Vec<BggPublicKey<T, P, M>>, BggError> {
        let log_q = self.sampler.modulus_bits();
        let columns = 2 * log_q * packed_input_size;
        let all_matrix = self
            .sampler
            .sample_hash(tag, 2, columns)
            .map_err(|e| BggError::SampleError(e.to_string()))?;
        let public_keys = (0..packed_input_size)
            .map(|idx| {
                let matrix = self
                    .matrix_op
                    .slice_rows(&all_matrix, 2 * log_q * idx, 2 * log_q * (idx + 1))
                    .expect("failed to slice the matrix");
                BggPublicKey { matrix, index: idx }
            })
            .collect_vec();
        Ok(public_keys)
    }
}

/// A sampler of an encoding in the BGG+ RLWE encoding scheme
///
/// # Fields
/// * `secret`: The secret vector.
/// * `error_sampler`: The sampler to generate LWE errors.
#[derive(Debug, Clone)]
pub struct BGGEncodingSampler<
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    S: PolyUniformSampler<T, P, M, GaussianDist>,
    G: PolyGadgetOps<T, P, M>,
> {
    pub(crate) secret_vec: PolyMatrix<T, P, M>,
    pub error_sampler: Arc<S>,
    pub poly_op: Arc<P>,
    pub matrix_op: Arc<M>,
    pub gadget_op: Arc<G>,
    _t: PhantomData<T>,
}

impl<
        T: PolyElemOps,
        P: PolyOps<T>,
        M: PolyMatrixOps<T, P>,
        S: PolyUniformSampler<T, P, M, GaussianDist>,
        G: PolyGadgetOps<T, P, M>,
    > BGGEncodingSampler<T, P, M, S, G>
{
    /// Create a new encoding sampler
    /// # Arguments
    /// * `secret`: The secret polynomial
    /// * `error_sampler`: The sampler to generate LWE errors
    /// * `matrix_op`: The matrix operations
    /// * `gadget_op`: The gadget operations
    pub fn new(
        secret: Poly<T, P>,
        error_sampler: Arc<S>,
        poly_op: Arc<P>,
        matrix_op: Arc<M>,
        gadget_op: Arc<G>,
    ) -> Self {
        let minus_one_poly = poly_op.minus_one().expect("failed to create a minus one polynomial");
        // 2*1 column vector
        let secret_vec = matrix_op.poly_vec_to_matrix(vec![secret.clone(), minus_one_poly]);
        // 1*2 row vector
        let secret_vec =
            matrix_op.transpose(&secret_vec).expect("failed to transpose the secret_vec");
        Self { secret_vec, error_sampler, poly_op, matrix_op, gadget_op, _t: PhantomData }
    }

    pub fn sample<R: RngCore>(
        &self,
        rng: &mut R,
        public_keys: &[BggPublicKey<T, P, M>],
        plaintexts: &[Poly<T, P>],
        reveal_plaintexts: bool,
    ) -> Result<Vec<BggEncoding<T, P, M>>, BggError> {
        let log_q = self.poly_op.modulus_bits();
        let packed_input_size = plaintexts.len();
        let columns = 2 * log_q * packed_input_size;
        let error = self
            .error_sampler
            .sample_uniform(rng, 1, columns)
            .map_err(|e| BggError::SampleError(e.to_string()))?;
        // [TODO] Avoid memory cloning here.
        let all_public_key_matrix = self
            .matrix_op
            .concat_rows(&public_keys.iter().map(|pk| pk.matrix.clone()).collect_vec())
            .map_err(|e| BggError::MatrixError(e.to_string()))?;
        let first_term = self
            .matrix_op
            .mul(&self.secret_vec, &all_public_key_matrix)
            .map_err(|e| BggError::MatrixError(e.to_string()))?;
        let gadget = self.gadget_op.gadget_matrix(2);
        let sg = self
            .matrix_op
            .mul(&self.secret_vec, &gadget)
            .map_err(|e| BggError::MatrixError(e.to_string()))?;
        let encoded_polys_vec = self.matrix_op.poly_vec_to_matrix(plaintexts.to_vec());
        let second_term = self
            .matrix_op
            .tensor(&encoded_polys_vec, &sg)
            .map_err(|e| BggError::MatrixError(e.to_string()))?;
        let all_vector = {
            let add1 = self
                .matrix_op
                .add(&first_term, &second_term)
                .map_err(|e| BggError::MatrixError(e.to_string()))?;
            self.matrix_op.add(&add1, &error).map_err(|e| BggError::MatrixError(e.to_string()))?
        };
        let encodings = plaintexts
            .iter()
            .enumerate()
            .map(|(idx, poly)| BggEncoding {
                vector: self
                    .matrix_op
                    .slice_columns(&all_vector, 2 * log_q * idx, 2 * log_q * (idx + 1))
                    .expect("failed to slice all_vector to sample the BGG+ encodings"),
                plaintext: if reveal_plaintexts { Some(poly.clone()) } else { None },
                index: idx,
            })
            .collect_vec();
        Ok(encodings)
    }
}
