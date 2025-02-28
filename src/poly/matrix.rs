use super::Poly;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg},
};

pub trait PolyMatrix:
    Sized + Clone + Debug + PartialEq + Eq + Add + Mul + Neg + Mul<Self::P>
{
    // type Error: std::error::Error + Send + Sync + 'static;
    type P: Poly;

    fn from_poly_vec(params: &<Self::P as Poly>::Params, vec: Vec<Vec<Self::P>>) -> Self;
    fn entry(&self, i: usize, j: usize) -> &Self::P;
    fn size(&self) -> (usize, usize);
    fn row_size(&self) -> usize {
        self.size().0
    }
    fn col_size(&self) -> usize {
        self.size().1
    }
    fn slice(
        &self,
        row_start: usize,
        row_end: usize,
        column_start: usize,
        column_end: usize,
    ) -> Self;
    fn slice_rows(&self, start: usize, end: usize) -> Self {
        let (_, columns) = self.size();
        self.slice(start, end, 0, columns)
    }
    fn slice_columns(&self, start: usize, end: usize) -> Self {
        let (rows, _) = self.size();
        self.slice(0, rows, start, end)
    }
    fn zero(params: &<Self::P as Poly>::Params, rows: usize, columns: usize) -> Self;
    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self;
    fn transpose(&self) -> Self;
    // (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[Self]) -> Self;
    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[Self]) -> Self;
    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[Self]) -> Self;
    fn tensor(&self, other: &Self) -> Self;
}
