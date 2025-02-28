use super::{Obfuscation, ObfuscationError};
use crate::bgg::*;
use crate::bgg::{eval::*, sampler::*, *};
use crate::poly::gadget::PolyGadgetOps;
use crate::poly::{matrix::*, sampler::*, *};
use crate::utils::*;
use num_traits::{One, Zero};
use rand::Rng;
use rand::RngCore;
use std::sync::Arc;

const TAG_R_0: &[u8] = b"R_0";
const TAG_R_1: &[u8] = b"R_1";
const TAG_A_FHE_BAR: &[u8] = b"A_FHE_BAR";
const TAG_BGG_PUBKEY_INPUT: &[u8] = b"BGG_PUBKEY_INPUT";
const TAG_BGG_PUBKEY_FHEKEY: &[u8] = b"BGG_PUBKEY_FHEKY";

#[derive(Debug, Clone)]
pub struct PublicSampledData<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    pub r_0: PolyMatrix<T, P, M>,
    pub r_1: PolyMatrix<T, P, M>,
    pub a_fhe_bar: PolyMatrix<T, P, M>,
    pub pubkeys_input: Vec<BggPublicKey<T, P, M>>,
    pub pubkeys_fhe_key: Vec<BggPublicKey<T, P, M>>,
    pub t_0: (PolyMatrix<T, P, M>, PolyMatrix<T, P, M>),
    pub t_1: (PolyMatrix<T, P, M>, PolyMatrix<T, P, M>),
}

impl<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> PublicSampledData<T, P, M> {
    pub fn sample<G: PolyGadgetOps<T, P, M>, S: PolyHashSampler<T, P, M>>(
        matrix_op: Arc<M>,
        gadget_op: Arc<G>,
        sampler: Arc<S>,
        packed_input_size: usize,
    ) -> Result<Self, ObfuscationError> {
        let r_0_bar = sampler
            .sample_hash::<_, BitDist>(TAG_R_0, 1, 1)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        let r_1_bar = sampler
            .sample_hash::<_, BitDist>(TAG_R_1, 1, 1)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        let one = matrix_op.identity(1, None);
        let r_0 = matrix_op
            .concat_diag(&[r_0_bar, one.clone()])
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let r_1 = matrix_op
            .concat_diag(&[r_1_bar, one.clone()])
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;

        let log_q = sampler.modulus_bits();
        let a_fhe_bar = sampler
            .sample_hash::<_, FinRingDist>(TAG_A_FHE_BAR, 2, 2 * log_q)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        let bgg_pubkey_sampler = BGGPublicKeySampler::new(sampler, matrix_op.clone());
        let pubkeys_input =
            bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_INPUT, packed_input_size + 1)?;
        let pubkeys_fhe_key = bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_FHEKEY, 2)?;
        let identity_input = matrix_op.identity(packed_input_size + 1, None);
        let gadget_2 = gadget_op.gadget_matrix(2);
        let identity_2 = matrix_op.identity(2, None);
        let mut ts = vec![];
        for bit in 0..1 {
            let r = if bit == 0 { r_0.clone() } else { r_1.clone() };
            let rg = matrix_op
                .mul(&r, &gadget_2)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let rg_decomposed = gadget_op
                .decompose(&rg)
                .map_err(|e| ObfuscationError::GadgetError(e.to_string()))?;
            let t_input = matrix_op
                .tensor(&identity_input, &rg_decomposed)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let t_fhe_key = matrix_op
                .tensor(&identity_2, &rg_decomposed)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            ts.push((t_input, t_fhe_key));
        }
        Ok(Self {
            r_0,
            r_1,
            a_fhe_bar,
            pubkeys_input,
            pubkeys_fhe_key,
            t_0: ts[0].clone(),
            t_1: ts[1].clone(),
        })
    }
}
