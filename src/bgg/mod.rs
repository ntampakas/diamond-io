// pub mod ciphertext;
// pub mod circuit;
// pub mod operations;
// pub mod parameters;
// pub mod eval;
// pub mod sampler;

use crate::poly::{matrix::*, *};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum BggError {
    #[error("Sample error: {0}")]
    SampleError(String),
    #[error("Poly error: {0}")]
    PolyError(String),
    #[error("Matrix error: {0}")]
    MatrixError(String),
    #[error("Gadget error: {0}")]
    GadgetError(String),
    #[error("Unknown plaintext for the left-hand input of multiplication: left-input index {0}, output index {1}")]
    UnknownPlaintextForMul(usize, usize),
}

#[derive(Debug, Clone)]
pub struct BggPublicKey<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    pub matrix: PolyMatrix<T, P, M>,
    pub index: usize,
}

#[derive(Debug, Clone)]
pub struct BggEncoding<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    pub vector: PolyMatrix<T, P, M>,
    pub plaintext: Option<Poly<T, P>>,
    pub index: usize,
}
