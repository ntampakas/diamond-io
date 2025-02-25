use openfhe::{cxx::CxxVector, ffi};

use super::{PolyDegree, PolyElemOps, PolyOps};

pub struct DCRTPoly<T>
where
    T: PolyElemOps,
{
    _phantom: Box<T>,
}

impl<T> PolyOps<T> for DCRTPoly<T>
where
    Self: PolyDegree<T>,
    T: PolyElemOps,
{
    type Error;

    type Poly;

    fn coeffs(&self, poly: &Self::Poly) -> &[super::PElem<T>] {
        todo!()
    }

    fn from_coeffs(coeffs: &[super::PElem<T>]) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn from_const(constant: &T) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn zero(&self) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn one(&self) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn minus_one(&self) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn neg(&self, a: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        todo!()
    }

    fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error> {
        todo!()
    }

    fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let minus_b = self.neg(b)?;
        self.add(a, &minus_b)
    }
}
