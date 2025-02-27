use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg},
};

use crate::poly::{Poly, PolyParams, PolynomialMatrix};

pub struct DCRTPolyMatrix<P>
where
    P: Poly,
{
    inner: Vec<Vec<P>>,
    params: P::Params,
    nrow: usize,
    ncol: usize,
}

// Add getter methods for inner and params
impl<P> DCRTPolyMatrix<P>
where
    P: Poly,
{
    pub fn inner(&self) -> &Vec<Vec<P>> {
        &self.inner
    }
    pub fn params(&self) -> &P::Params {
        &self.params
    }
}

impl<P> PolynomialMatrix<P> for DCRTPolyMatrix<P>
where
    P: Poly<Params = PolyParams> + Mul<Output = P> + Neg<Output = P> + 'static,
{
    type Error = std::io::Error;
    type Matrix = Vec<Vec<P>>;

    fn from_poly_vec(params: &P::Params, vec: Vec<Vec<P>>) -> Self {
        let mut c = vec![vec![P::const_zero(params); vec[0].len()]; vec.len()];
        for (i, row) in vec.iter().enumerate() {
            for (j, element) in row.iter().enumerate() {
                c[i][j] = element.clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow: vec.len(), ncol: vec[0].len() }
    }

    fn row_size(&self) -> usize {
        self.nrow
    }

    fn col_size(&self) -> usize {
        self.ncol
    }

    fn identity(&self) -> Result<Self, Self::Error> {
        let nrow = self.row_size();
        let ncol = self.col_size();
        if nrow != ncol {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Identity matrix must be square (ROWS must equal COLUMNS)",
            ));
        }
        let mut result = [vec![P::const_zero(&self.params)]];
        for i in 0..nrow {
            result[i][i] = P::const_one(&self.params);
        }

        Ok(Self { inner: result.to_vec(), params: self.params.clone(), ncol, nrow })
    }

    fn zero(params: &P::Params, nrow: usize, ncol: usize) -> Self {
        let mut c: Vec<Vec<P>> = vec![vec![P::const_zero(params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                c[i][j] = P::const_zero(params).clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow, ncol }
    }
}

// ====== Arithmetic ======

impl<P> Add for DCRTPolyMatrix<P>
where
    P: Poly<Params = PolyParams> + Mul<Output = P> + Neg<Output = P> + 'static,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let nrow = self.row_size();
        let ncol = self.col_size();
        let mut result = self.inner;
        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] += rhs.inner[i][j].clone();
            }
        }
        Self { inner: result.to_vec(), params: self.params, ncol, nrow }
    }
}

impl<P> Neg for DCRTPolyMatrix<P>
where
    P: Poly + Neg<Output = P> + 'static,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        let mut c: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); self.ncol]; self.nrow];
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                c[i][j] = -self.inner[i][j].clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: self.params, nrow: self.nrow, ncol: self.ncol }
    }
}

impl<P> Mul for DCRTPolyMatrix<P>
where
    P: Poly + Mul<Output = P> + 'static,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let nrow = self.nrow;
        let ncol = rhs.ncol;
        if rhs.nrow != self.ncol {
            panic!(
                "Multiplication condition failed: rhs.nrow ({}) must equal self.ncol ({})",
                rhs.nrow, self.ncol
            );
        }
        let common = self.ncol;
        let mut c: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                for k in 0..common {
                    c[i][j] += self.inner[i][k].clone() * rhs.inner[k][j].clone();
                }
            }
        }
        DCRTPolyMatrix { inner: c, params: self.params, nrow, ncol }
    }
}

// ====== Traits ======

impl<P> Clone for DCRTPolyMatrix<P>
where
    P: Poly<Params = PolyParams>,
{
    fn clone(&self) -> Self {
        Self {
            inner: self.inner.clone(),
            params: self.params.clone(),
            nrow: self.nrow,
            ncol: self.ncol,
        }
    }
}

impl<P> Debug for DCRTPolyMatrix<P>
where
    P: Poly<Params = PolyParams>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix").field("inner", &self.inner).finish()
    }
}

impl<P> PartialEq for DCRTPolyMatrix<P>
where
    P: Poly<Params = PolyParams>,
{
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<P> Eq for DCRTPolyMatrix<P> where P: Poly<Params = PolyParams> {}
