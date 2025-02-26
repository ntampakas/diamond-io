use super::Polynomial;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg},
};

pub trait PolynomialMatrix<P: Polynomial>:
    Sized + Clone + Debug + PartialEq + Eq + Add + Mul + Neg
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Matrix;

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
    fn zero(params: &P::Params, rows: usize, columns: usize) -> Self;
    fn from_poly_vec(params: &P::Params, vec: Self::Matrix) -> Self;
    // fn identity(&self, size: usize, scale: Option<&P>) -> Result<Self::Matrix, Self::Error>;
    fn identity(&self) -> Result<Self, Self::Error>;
    // fn transpose(&self, matrix: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
    // // (m * n1), (m * n2) -> (m * (n1 + n2))
    // fn concat_columns(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    // fn concat_rows(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    // fn concat_diag(&self, matrices: &[Self::Matrix]) -> Result<Self::Matrix, Self::Error>;
    // fn scalar_mul(&self, a: &Self::Matrix, scalar: &P) -> Result<Self::Matrix, Self::Error>;
    // fn tensor(&self, a: &Self::Matrix, b: &Self::Matrix) -> Result<Self::Matrix, Self::Error>;
}
