use std::{
    fmt::Debug,
    ops::{Add, Neg},
};

use nalgebra::SMatrix;
use num_traits::Zero;

use super::{matrix::PolynomialMatrix, Polynomial};

pub struct DCRTPolyMatrix<P, const ROWS: usize, const COLUMNS: usize>
where
    P: Polynomial,
{
    pub inner: SMatrix<P, ROWS, COLUMNS>,
}

impl<P, const ROWS: usize, const COLUMNS: usize> PolynomialMatrix<P>
    for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    type Error = std::io::Error;

    fn from_slice(slice: &[P]) -> Self {
        let mut z = SMatrix::<P, ROWS, COLUMNS>::zeros();
        z.copy_from_slice(slice);
        Self { inner: z }
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Clone for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    fn clone(&self) -> Self {
        Self { inner: self.inner.clone() }
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

impl<P, const ROWS: usize, const COLUMNS: usize> Add for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    type Output = Self;
    fn add(self, _rhs: Self) -> Self {
        todo!()
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Neg for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        todo!()
    }
}

impl<P, const ROWS: usize, const COLUMNS: usize> Zero for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: Polynomial + 'static,
{
    fn zero() -> Self {
        Self { inner: SMatrix::<P, ROWS, COLUMNS>::zeros() }
    }

    fn is_zero(&self) -> bool {
        todo!()
    }
}
