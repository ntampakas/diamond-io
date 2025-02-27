use super::Poly;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg},
};

pub trait PolyMatrix: Sized + Clone + Debug + PartialEq + Eq + Add + Mul + Neg {
    type Error: std::error::Error + Send + Sync + 'static;
    type P: Poly;

    fn from_poly_vec(
        params: &Self::P::Params,
        vec: &[P],
        rows: usize,
        columns: usize,
    ) -> Result<Self, Self::Error>;
    fn entry(&self, i: usize, j: usize) -> P;
    fn size(&self) -> (usize, usize);
    fn row_size(&self) -> usize {
        self.size().0
    }
    fn col_size(&self) -> usize {
        self.size().1
    }
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
    fn zero(params: &Self::P::Params, rows: usize, columns: usize) -> Self;
    fn identity(params: &Self::P::Params, size: usize, scalar: Option<P>) -> Self;
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
