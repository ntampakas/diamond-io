use crate::poly::{Poly, PolyMatrix, PolyParams};
use openfhe::ffi::DCRTPolyParams;
use std::{
    fmt::Debug,
    ops::{Add, Mul, Neg},
};
// use thiserror::Error;

// #[derive(Debug, Error)]
// pub enum DCRTPolyMatrixError {
//     #[error("Concat error: while the shape of the first matrix is ({0}, {1}), that of the {2}-th matirx is ({3},{4})")]
//     ConcatError(usize, usize, usize, usize, usize),
// }

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

impl<P> PolyMatrix for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Neg<Output = P> + 'static,
{
    // type Error = DCRTPolyMatrixError;
    type P = P;

    fn from_poly_vec(params: &P::Params, vec: Vec<Vec<P>>) -> Self {
        let mut c = vec![vec![P::const_zero(params); vec[0].len()]; vec.len()];
        for (i, row) in vec.iter().enumerate() {
            for (j, element) in row.iter().enumerate() {
                c[i][j] = element.clone();
            }
        }
        DCRTPolyMatrix { inner: c, params: params.clone(), nrow: vec.len(), ncol: vec[0].len() }
    }

    fn entry(&self, i: usize, j: usize) -> &Self::P {
        &self.inner[i][j]
    }

    fn size(&self) -> (usize, usize) {
        (self.nrow, self.ncol)
    }

    fn row_size(&self) -> usize {
        self.nrow
    }

    fn col_size(&self) -> usize {
        self.ncol
    }

    fn slice(
        &self,
        row_start: usize,
        row_end: usize,
        column_start: usize,
        column_end: usize,
    ) -> Self {
        let mut c =
            vec![vec![P::const_zero(&self.params); column_end - column_start]; row_end - row_start];
        for i in row_start..row_end {
            for j in column_start..column_end {
                c[i - row_start][j - column_start] = self.inner[i][j].clone();
            }
        }
        DCRTPolyMatrix {
            inner: c,
            params: self.params.clone(),
            nrow: row_end - row_start,
            ncol: column_end - column_start,
        }
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

    fn identity(params: &<Self::P as Poly>::Params, size: usize, scalar: Option<Self::P>) -> Self {
        let nrow = size;
        let ncol = size;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(params); ncol]; nrow];
        let scalr = scalar.unwrap_or(P::const_one(params));
        for i in 0..size {
            result[i][i] = scalar;
        }
        DCRTPolyMatrix { inner: result, params: params.clone(), nrow, ncol }
    }

    fn transpose(&self) -> Self {
        let nrow = self.ncol;
        let ncol = self.nrow;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];
        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[j][i].clone();
            }
        }
        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    // (m * n1), (m * n2) -> (m * (n1 + n2))
    fn concat_columns(&self, others: &[Self]) -> Self {
        for (idx, other) in others.iter().enumerate() {
            if self.nrow != other.nrow {
                panic!("Concat error: while the shape of the first matrix is ({0}, {1}), that of the {2}-th matirx is ({3},{4})",self.nrow,self.ncol,idx,other.nrow,other.ncol);
            }
        }
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; self.nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut offset = self.ncol;
        for other in others {
            for i in 0..self.nrow {
                for j in 0..other.ncol {
                    result[i][offset + j] = other.inner[i][j].clone();
                }
            }
            offset += other.ncol;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow: self.nrow, ncol }
    }

    // (m1 * n), (m2 * n) -> ((m1 + m2) * n)
    fn concat_rows(&self, others: &[Self]) -> Self {
        for (idx, other) in others.iter().enumerate() {
            if self.ncol != other.ncol {
                panic!("Concat error: while the shape of the first matrix is ({0}, {1}), that of the {2}-th matirx is ({3},{4})",self.nrow,self.ncol,idx,other.nrow,other.ncol);
            }
        }
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); self.ncol]; nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut offset = self.nrow;
        for other in others {
            for i in 0..other.nrow {
                for j in 0..other.ncol {
                    result[offset + i][j] = other.inner[i][j].clone();
                }
            }
            offset += other.nrow;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol: self.ncol }
    }

    // (m1 * n1), (m2 * n2) -> ((m1 + m2) * (n1 + n2))
    fn concat_diag(&self, others: &[Self]) -> Self {
        let nrow = self.nrow + others.iter().map(|x| x.nrow).sum::<usize>();
        let ncol = self.ncol + others.iter().map(|x| x.ncol).sum::<usize>();
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        // Copy elements from self
        for i in 0..self.nrow {
            for j in 0..self.ncol {
                result[i][j] = self.inner[i][j].clone();
            }
        }

        // Copy elements from others
        let mut row_offset = self.nrow;
        let mut col_offset = self.ncol;
        for other in others {
            for i in 0..other.nrow {
                for j in 0..other.ncol {
                    result[row_offset + i][col_offset + j] = other.inner[i][j].clone();
                }
            }
            row_offset += other.nrow;
            col_offset += other.ncol;
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }

    fn tensor(&self, other: &Self) -> Self {
        let nrow = self.nrow * other.nrow;
        let ncol = self.ncol * other.ncol;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        for i1 in 0..self.nrow {
            for j1 in 0..self.ncol {
                for i2 in 0..other.nrow {
                    for j2 in 0..other.ncol {
                        let i = i1 * other.nrow + i2;
                        let j = j1 * other.ncol + j2;
                        result[i][j] = self.inner[i1][j1].clone() * other.inner[i2][j2].clone();
                    }
                }
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// ====== Arithmetic ======

impl<P> Add for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Neg<Output = P> + 'static,
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
    P: Poly<Params = DCRTPolyParams> + Neg<Output = P> + 'static,
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
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + 'static,
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

// Implement multiplication of a matrix by a polynomial
impl<P> Mul<P> for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Clone + 'static,
{
    type Output = Self;

    fn mul(self, rhs: P) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[i][j].clone() * rhs.clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params, nrow, ncol }
    }
}

// Implement multiplication of a matrix reference by a polynomial
impl<P> Mul<P> for &DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Clone + 'static,
{
    type Output = DCRTPolyMatrix<P>;

    fn mul(self, rhs: P) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[i][j].clone() * rhs.clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// Implement multiplication of a matrix by a polynomial reference
impl<P> Mul<&P> for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Clone + 'static,
{
    type Output = Self;

    fn mul(self, rhs: &P) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[i][j].clone() * rhs.clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params, nrow, ncol }
    }
}

// Implement multiplication of a matrix reference by a polynomial reference
impl<P> Mul<&P> for &DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams> + Mul<Output = P> + Clone + 'static,
{
    type Output = DCRTPolyMatrix<P>;

    fn mul(self, rhs: &P) -> Self::Output {
        let nrow = self.nrow;
        let ncol = self.ncol;
        let mut result: Vec<Vec<P>> = vec![vec![P::const_zero(&self.params); ncol]; nrow];

        for i in 0..nrow {
            for j in 0..ncol {
                result[i][j] = self.inner[i][j].clone() * rhs.clone();
            }
        }

        DCRTPolyMatrix { inner: result, params: self.params.clone(), nrow, ncol }
    }
}

// ====== Traits ======

impl<P> Clone for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams>,
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
    P: Poly<Params = DCRTPolyParams>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DCRTPolyMatrix").field("inner", &self.inner).finish()
    }
}

impl<P> PartialEq for DCRTPolyMatrix<P>
where
    P: Poly<Params = DCRTPolyParams>,
{
    fn eq(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}

impl<P> Eq for DCRTPolyMatrix<P> where P: Poly<Params = DCRTPolyParams> {}
