use std::{
    fmt::Debug,
    ops::{Add, Neg},
};

use num_traits::Zero;

use super::{
    matrix::{get_null_matrix, get_zero_matrix, Matrix, PolynomialMatrix},
    params::Params,
    Polynomial,
};

/// matrix multiplication
#[allow(clippy::needless_range_loop)]
pub fn mult<P, const COMMON: usize, const R1: usize, const R2: usize>(
    a: &DCRTPolyMatrix<P, R1, COMMON>,
    b: &DCRTPolyMatrix<P, COMMON, R2>,
    params: Params,
) -> DCRTPolyMatrix<P, R1, R2>
where
    P: Polynomial + 'static,
{
    let mut c = get_zero_matrix::<P, R1, R2>(&params);
    for i in 0..R1 {
        for j in 0..R2 {
            for k in 0..COMMON {
                c[i][j] += a.inner[i][k].clone() * b.inner[k][j].clone();
            }
        }
    }
    DCRTPolyMatrix { inner: c, params: Some(params) }
}

pub struct DCRTPolyMatrix<P, const ROW: usize, const COLUMNS: usize>
where
    P: Polynomial,
{
    pub inner: Matrix<P, ROW, COLUMNS>,
    pub params: Option<Params>,
}

impl<P, const ROWS: usize, const COLUMNS: usize> PolynomialMatrix<P, ROWS, COLUMNS>
    for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    type Error = std::io::Error;

    fn from_slice(slice: &[[P; COLUMNS]; ROWS]) -> Self {
        let mut c = get_null_matrix::<P, ROWS, COLUMNS>();
        for (i, row) in slice.iter().enumerate() {
            for (j, element) in row.iter().enumerate() {
                c[i][j] = element.clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: None }
    }

    fn row_size(&self) -> usize {
        ROWS
    }

    fn col_size(&self) -> usize {
        COLUMNS
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Clone for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    fn clone(&self) -> Self {
        Self { inner: self.inner.clone(), params: self.params.clone() }
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Debug for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix").field("inner", &self.inner).finish()
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> PartialEq for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Eq for DCRTPolyMatrix<P, ROWS, COLUMNS> where
    P: Polynomial
{
}

#[allow(clippy::needless_range_loop)]
impl<P, const ROWS: usize, const COLUMNS: usize> Add for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = self.inner;
        for i in 0..ROWS {
            for j in 0..COLUMNS {
                result[i][j] += rhs.inner[i][j].clone();
            }
        }
        Self { inner: result, params: self.params }
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Neg for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        todo!()
        // let mut c = get_zero_matrix::<P, ROWS, COLUMNS>(&self.params.unwrap());
        // for i in 0..ROWS {
        //     for j in 0..COLUMNS {
        //         let neg: P = -self.inner[i][j];
        //         c[i][j] = neg;
        //     }
        // }
        // Self { inner: c, params: self.params }
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Zero for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    fn zero() -> Self {
        todo!()
    }

    fn is_zero(&self) -> bool {
        todo!()
    }
}

// impl<P, const ROWS: usize, const COLUMNS: usize> One for DCRTPolyMatrix<P, ROWS, COLUMNS>
// where
//     P: Polynomial + 'static,
// {
//     fn one() -> Self {
//         todo!()
//     }
// }
