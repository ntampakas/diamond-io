use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPolyImpl},
};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub},
    sync::Arc,
};

use crate::poly::{PElem, PolyParams, Polynomial};

// ====== FieldElement ======

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, PartialOrd, Ord)]
pub struct FieldElement {
    value: u64,   // TODO: support BigInt
    modulus: u64, // TODO: support BigInt
}

impl FieldElement {
    pub fn new(value: u64, modulus: u64) -> Self {
        let reduced_value = value % modulus;
        Self { value: reduced_value, modulus }
    }

    pub fn value(&self) -> u64 {
        self.value
    }

    pub fn modulus(&self) -> u64 {
        self.modulus
    }
}

// ====== FieldElement ======

impl Add for FieldElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new((self.value + rhs.value) % self.modulus, self.modulus)
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new((self.value * rhs.value) % self.modulus, self.modulus)
    }
}

impl PElem for FieldElement {}

// ====== DCRTPoly ======

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
