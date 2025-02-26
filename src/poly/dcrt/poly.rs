use openfhe::{
    cxx::{CxxVector, UniquePtr},
    ffi::{self, DCRTPolyImpl},
};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub},
    sync::Arc,
};

use crate::poly::{PolyParams, Polynomial};

use super::FieldElement;

#[derive(Clone, Debug)]
pub struct DCRTPoly {
    pub ptr_poly: Arc<UniquePtr<DCRTPolyImpl>>,
}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyImpl>) -> Self {
        Self { ptr_poly: ptr_poly.into() }
    }
}

impl Polynomial for DCRTPoly {
    type Error = std::io::Error;
    type Elem = FieldElement;
    type Params = PolyParams;

    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Result<Self, Self::Error> {
        // TODO: check if coeffs modulus is the same as params modulus
        // TODO: check if coeffs length is the same as the ring size
        let mut coeffs_cxx = CxxVector::<i64>::new();
        for coeff in coeffs {
            coeffs_cxx.pin_mut().push(coeff.value().try_into().unwrap());
        }
        let res = ffi::DCRTPolyGenFromVec(&params.ptr_params, &coeffs_cxx);
        Ok(DCRTPoly::new(res))
    }

    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Result<Self, Self::Error> {
        let res = ffi::DCRTPolyGenFromConst(&params.ptr_params, constant.value());
        Ok(DCRTPoly::new(res))
    }

    fn const_zero(params: &Self::Params) -> Self {
        let res = ffi::DCRTPolyGenFromConst(&params.ptr_params, 0);
        DCRTPoly::new(res)
    }

    fn const_one(params: &Self::Params) -> Self {
        let res = ffi::DCRTPolyGenFromConst(&params.ptr_params, 1);
        DCRTPoly::new(res)
    }
}

// ====== Arithmetic ======

impl Add for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly);
        DCRTPoly::new(res)
    }
}

impl Mul for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly);
        DCRTPoly::new(res)
    }
}

impl Sub for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        let minus_rhs = rhs.neg();
        self.add(minus_rhs)
    }
}

impl Neg for DCRTPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let res = self.ptr_poly.Negate();
        DCRTPoly::new(res)
    }
}

impl PartialEq for DCRTPoly {
    fn eq(&self, other: &Self) -> bool {
        if self.ptr_poly.is_null() || other.ptr_poly.is_null() {
            return self.ptr_poly.is_null() && other.ptr_poly.is_null();
        }
        self.ptr_poly.IsEqual(&other.ptr_poly)
    }
}

impl Eq for DCRTPoly {}

impl AddAssign for DCRTPoly {
    fn add_assign(&mut self, rhs: Self) {
        if self.ptr_poly.is_null() || rhs.ptr_poly.is_null() {
            panic!("Attempted to dereference a null pointer");
        }
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl MulAssign for DCRTPoly {
    fn mul_assign(&mut self, rhs: Self) {
        if self.ptr_poly.is_null() || rhs.ptr_poly.is_null() {
            panic!("Attempted to dereference a null pointer");
        }
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}
