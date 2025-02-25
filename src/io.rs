pub mod eval;
pub mod obf;
pub mod utils;

use crate::bgg::{BggEncoding, BggError, BggPublicKey};
use crate::poly::gadget::PolyGadgetOps;
use crate::poly::{matrix::*, sampler::*, *};
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
{
    pub hash_key: Vec<u8>,
    // pub fhe_enc: PolyMatrix<T, P, M>,
    pub b_fhe: PolyMatrix<T, P, M>,
    pub encode_input: Vec<BggEncoding<T, P, M>>,
    pub encode_fhe_key: Vec<BggEncoding<T, P, M>>,
    pub p_init: PolyMatrix<T, P, M>,
    pub m_preimages: Vec<(PolyMatrix<T, P, M>, PolyMatrix<T, P, M>)>,
    pub n_preimages: Vec<(PolyMatrix<T, P, M>, PolyMatrix<T, P, M>)>,
    pub k_preimages: Vec<(PolyMatrix<T, P, M>, PolyMatrix<T, P, M>)>,
    // pub final_preimage: PolyMatrix<T, P, M>,
}
