use num_traits::Zero;

use super::{params::Params, Polynomial};
use std::{
    array,
    fmt::Debug,
    ops::{Add, Neg},
};

pub type Matrix<P, const ROW: usize, const COL: usize> = [[P; COL]; ROW];

pub fn get_zero_matrix<P, const ROW: usize, const COL: usize>(
    params: &Params,
) -> Matrix<P, ROW, COL>
where
    P: Polynomial,
{
    array::from_fn(|_| array::from_fn(|_| P::const_zero(params)))
}

pub fn get_null_matrix<P, const ROW: usize, const COL: usize>() -> Matrix<P, ROW, COL>
where
    P: Polynomial,
{
    array::from_fn(|_| array::from_fn(|_| P::zero()))
}

// pub type PolynomialMatrix<P, M> = <M as PolyMatrixOps<P>>::Matrix;

pub trait PolynomialMatrix<P: Polynomial, const ROWS: usize, const COLUMNS: usize>:
    Sized + Clone + Debug + PartialEq + Eq + Add + Neg + Zero
{
    type Error: std::error::Error + Send + Sync + 'static;
    // fn poly_vec_to_matrix(&self, polys: Vec<P>) -> Self::Matrix;
    fn row_size(&self) -> usize;
    fn col_size(&self) -> usize;
    // fn entry(&self, matrix: &Self::Matrix, i: usize, j: usize) -> Result<P, Self::Error>;
    // fn size(&self, matrix: &Self::Matrix) -> Result<(usize, usize), Self::Error>;
    // fn slice(
    //     &self,
    //     matrix: &Self::Matrix,
    //     row_start: usize,
    //     row_end: usize,
    //     column_start: usize,
    //     column_end: usize,
    // ) -> Result<Self::Matrix, Self::Error>;
    // fn slice_rows(
    //     &self,
    //     matrix: &Self::Matrix,
    //     start: usize,
    //     end: usize,
    // ) -> Result<Self::Matrix, Self::Error> {
    //     let (_, columns) = self.size(matrix)?;
    //     self.slice(matrix, start, end, 0, columns)
    // }
    // fn slice_columns(
    //     &self,
    //     matrix: &Self::Matrix,
    //     start: usize,
    //     end: usize,
    // ) -> Result<Self::Matrix, Self::Error> {
    //     let (rows, _) = self.size(matrix)?;
    //     self.slice(matrix, 0, rows, start, end)
    // }
    // fn zero() -> Self;
    fn from_slice(slice: &[[P; COLUMNS]; ROWS]) -> Self;
    // fn identity(&self, size: usize, scale: Option<&P>) -> Result<Self::Matrix, Self::Error>;
    // fn transpose(&self, matrix: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
    // // (m * n1), (m * n2) -> (m * (n1 + n2))
    // fn concat_columns(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    // fn concat_rows(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    // fn concat_diag(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // fn neg(&self, a: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
    // fn sub(&self, a: &Self::Matrix, b: &Self::Matrix) -> Result<Self::Matrix, Self::Error> {
    //     let minus_b = self.neg(b)?;
    //     self.add(a, &minus_b)
    // }
    // fn scalar_mul(&self, a: &Self::Matrix, scalar: &P) -> Result<Self::Matrix, Self::Error>;
    // fn mul(&self, a: &Self, b: &Self) -> Result<Self, Self::Error>;
    // fn tensor(&self, a: &Self::Matrix, b: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
}
