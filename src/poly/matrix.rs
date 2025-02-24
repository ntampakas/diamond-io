use super::{Poly, PolyElemOps, PolyOps};

pub type PolyMatrix<T, P, M> = <M as PolyMatrixOps<T, P>>::Matrix;

pub trait PolyMatrixOps<T: PolyElemOps, P: PolyOps<T>> {
    type Error: std::error::Error;
    type Matrix;
    fn entry(&self, matrix: &Self::Matrix, i: usize, j: usize) -> Result<Poly<T, P>, Self::Error>;
    fn transpose(&self, matrix: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
    fn add(&self, a: &Self::Matrix, b: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
    fn mul(&self, a: &Self::Matrix, b: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
}
