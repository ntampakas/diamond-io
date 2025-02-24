use super::{
    matrix::{PolyMatrix, PolyMatrixOps},
    Poly, PolyElemOps, PolyOps,
};

pub trait PolyGadgetOps<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    type Error: std::error::Error;
    fn gadget_vector(&self) -> PolyMatrix<T, P, M>;
    fn gadget_matrix(&self, n: usize) -> PolyMatrix<T, P, M>;
    fn decompose(&self, matrix: &PolyMatrix<T, P, M>) -> Result<PolyMatrix<T, P, M>, Self::Error>;
}
