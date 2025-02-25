use nalgebra::SMatrix;
use std::fmt::Debug;

use super::{matrix::PolyMatrixOps, PolyOps};

pub struct DCRTPolyMatrix<P, const ROWS: usize, const COLUMNS: usize>
where
    P: PolyOps,
{
    pub inner: SMatrix<P, ROWS, COLUMNS>,
}

impl<P, const ROWS: usize, const COLUMNS: usize> PolyMatrixOps<P>
    for DCRTPolyMatrix<P, ROWS, COLUMNS>
where
    P: PolyOps + Clone + PartialEq + Debug + num_traits::identities::Zero + 'static,
{
    type Error = std::io::Error;
    type Matrix = Self;

    fn add(&self, _rhs: &Self::Matrix, _lhs: &Self::Matrix) -> Result<Self::Matrix, Self::Error> {
        todo!()
    }

    fn mul(&self, _rhs: &Self::Matrix, _lhs: &Self::Matrix) -> Result<Self::Matrix, Self::Error> {
        todo!()
    }

    fn zero() -> Self::Matrix {
        Self { inner: SMatrix::<P, ROWS, COLUMNS>::zeros() }
    }

    fn from_slice(slice: &[P]) -> Self::Matrix {
        let mut z = SMatrix::<P, ROWS, COLUMNS>::zeros();
        z.copy_from_slice(slice);
        Self { inner: z }
    }
}
