// pub mod eval;
// pub mod obf;

use crate::{
    bgg::{BggEncoding, BggError, BggPublicKey},
    poly::{matrix::*, *},
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum ObfuscationError {
    #[error("Sample error: {0}")]
    SampleError(String),
    #[error("Poly error: {0}")]
    PolyError(String),
    #[error("Matrix error: {0}")]
    MatrixError(String),
    #[error("Gadget error: {0}")]
    GadgetError(String),
    #[error(transparent)]
    BggError(#[from] BggError),
}

#[derive(Debug, Clone)]
pub struct Obfuscation<T, P, M>
where
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    // S: PolyUniformSampler<T, P, M, BitDist>
    //     + PolyUniformSampler<T, P, M, GaussianDist>
    //     + PolyUniformSampler<T, P, M, FinRingDist>
    //     + PolyUniformSampler<T, P, M, BitDist>,
{
    pub hash_key: Vec<u8>,
    pub fhe_enc: PolyMatrix<T, P, M>,
    pub input_encode: BggEncoding<T, P, M>,
    pub fhe_key_encode: BggPublicKey<T, P, M>,
    pub init_p: PolyMatrix<T, P, M>,
    pub m_preimages: PolyMatrix<T, P, M>,
    pub n_preimages: PolyMatrix<T, P, M>,
    pub k_preimages: PolyMatrix<T, P, M>,
    pub final_preimage: PolyMatrix<T, P, M>,
}
