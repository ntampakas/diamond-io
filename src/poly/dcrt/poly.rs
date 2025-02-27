use openfhe::{
    cxx::{CxxVector, UniquePtr},
    ffi::{self, DCRTPolyImpl},
};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub},
    sync::Arc,
};

use super::fin_ring::FinRing;
use super::params::DCRTPolyParams;
use crate::poly::{Poly, PolyElem, PolyParams};

#[derive(Clone, Debug)]
pub struct DCRTPoly {
    ptr_poly: Arc<UniquePtr<DCRTPolyImpl>>,
}

impl DCRTPoly {
    pub fn new(ptr_poly: UniquePtr<DCRTPolyImpl>) -> Self {
        Self { ptr_poly: ptr_poly.into() }
    }

    pub fn get_poly(&self) -> &UniquePtr<DCRTPolyImpl> {
        &self.ptr_poly
    }
}

impl Poly for DCRTPoly {
    type Elem = FinRing;
    type Params = DCRTPolyParams;

    fn coeffs(&self) -> Vec<Self::Elem> {
        todo!()
    }

    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Self {
        // TODO: check if coeffs modulus is the same as params modulus
        // TODO: check if coeffs length is the same as the ring size
        let mut coeffs_cxx = CxxVector::<i64>::new();
        for coeff in coeffs {
            coeffs_cxx.pin_mut().push(coeff.value().try_into().unwrap());
        }
        let res = ffi::DCRTPolyGenFromVec(params.get_params(), &coeffs_cxx);
        DCRTPoly::new(res)
    }

    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self {
        let res =
            ffi::DCRTPolyGenFromConst(params.get_params(), constant.value().to_u64_digits()[0]);
        DCRTPoly::new(res)
    }

    fn const_zero(params: &Self::Params) -> Self {
        let res = ffi::DCRTPolyGenFromConst(params.get_params(), 0);
        DCRTPoly::new(res)
    }

    fn const_one(params: &Self::Params) -> Self {
        let res = ffi::DCRTPolyGenFromConst(params.get_params(), 1);
        DCRTPoly::new(res)
    }

    fn const_minus_one(params: &Self::Params) -> Self {
        // let fe = FinRing::minus_one(params);
        // let res = ffi::DCRTPolyGenFromConst(params.get_params(), fe.value().to_u64_digits()[0]);
        // Ok(DCRTPoly::new(res))
        todo!()
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

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_dcrtpoly_arithmetic() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let q = Arc::new(params.modulus());

        // todo: replace value and modulus from param
        let coeffs1 = [
            FinRing::new(100u32, q.clone()),
            FinRing::new(200u32, q.clone()),
            FinRing::new(300u32, q.clone()),
            FinRing::new(400u32, q.clone()),
        ];
        let coeffs2 = [
            FinRing::new(500u32, q.clone()),
            FinRing::new(600u32, q.clone()),
            FinRing::new(700u32, q.clone()),
            FinRing::new(800u32, q.clone()),
        ];

        // 3. Create polynomials from those coefficients.
        let poly1 = DCRTPoly::from_coeffs(&params, &coeffs1);
        let poly2 = DCRTPoly::from_coeffs(&params, &coeffs2);

        // 4. Test addition.
        let sum = poly1.clone() + poly2.clone();

        // 5. Test multiplication.
        let product = poly1.clone() * poly2.clone();

        // 6. Test negation / subtraction.
        let neg_poly2 = poly2.clone().neg();
        let difference = poly1.clone() - poly2.clone();

        let mut poly_add_assign = poly1.clone();
        poly_add_assign += poly2.clone();

        let mut poly_mul_assign = poly1.clone();
        poly_mul_assign *= poly2.clone();

        // 8. Make some assertions
        assert!(sum != poly1, "Sum should differ from original poly1");
        assert!(neg_poly2 != poly2, "Negated polynomial should differ from original");
        assert_eq!(difference + poly2, poly1, "p1 - p2 + p2 should be p1");

        assert_eq!(poly_add_assign, sum, "+= result should match separate +");
        assert_eq!(poly_mul_assign, product, "*= result should match separate *");

        // 9. Test from_const / const_zero / const_one
        // todo: `get_coeffs``
        // let const_poly = DCRTPoly::from_const(&params, &FieldElement::new(123, dummy_modulus))
        //     .expect("Failed to create DCRTPoly from const");
        // assert_eq!(
        //     const_poly,
        //     DCRTPoly::from_coeffs(&params, &[FieldElement::new(123, dummy_modulus); 1]).unwrap(),
        //     "from_const should produce a polynomial with all coeffs = 123"
        // );

        let zero_poly = DCRTPoly::const_zero(&params);
        assert_eq!(zero_poly, zero_poly.clone() + zero_poly.clone(), "0 + 0 = 0");

        let one_poly = DCRTPoly::const_one(&params);
        assert_eq!(zero_poly, one_poly.clone() - one_poly.clone(), "1 - 1 = 0");
    }
}
