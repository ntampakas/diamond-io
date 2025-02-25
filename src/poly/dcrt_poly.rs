use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPolyImpl},
};

use super::PolyOps;

pub struct DCRTPoly {
    ptr_poly: UniquePtr<DCRTPolyImpl>,
}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyImpl>) -> Self {
        Self { ptr_poly }
    }
}

impl PolyOps for DCRTPoly {
    type Error = std::io::Error;
    type Poly = Self;

    fn add(&self, rhs: &Self::Poly, lhs: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &lhs.ptr_poly);
        Ok(DCRTPoly::new(res))
    }

    fn mul(&self, rhs: &Self::Poly, lhs: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &lhs.ptr_poly);
        Ok(DCRTPoly::new(res))
    }
}
