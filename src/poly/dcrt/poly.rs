use super::{
    element::{FinRingElem, FinRingElemError},
    params::DCRTPolyParams,
};
use crate::poly::{Poly, PolyParams};
use num_bigint::BigUint;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, DCRTPolyImpl},
};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};
use thiserror::Error;

#[derive(Error, Debug)]
pub enum DCRTPolyError {
    #[error("fin ring error: {0}")]
    FinRingElem(#[from] FinRingElemError),

    #[error("unmatched modulus error: expected {expected}, found {found}")]
    UnmatchedModulus { expected: String, found: String },
}

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
    type Elem = FinRingElem;
    type Params = DCRTPolyParams;
    type Error = DCRTPolyError;

    fn coeffs(&self) -> Result<Vec<Self::Elem>, Self::Error> {
        let coeffs = self
            .ptr_poly
            .GetCoefficients()
            .iter()
            .map(|s| {
                FinRingElem::from_str(s, &self.ptr_poly.GetModulus())
                    .map_err(DCRTPolyError::FinRingElem)
            })
            .collect();
        coeffs
    }

    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Result<Self, Self::Error> {
        let mut coeffs_cxx = Vec::with_capacity(coeffs.len());
        let modulus = params.modulus();
        for coeff in coeffs {
            let coeff_modulus = coeff.modulus();
            if coeff_modulus != modulus.as_ref() {
                return Err(DCRTPolyError::UnmatchedModulus {
                    expected: modulus.to_string(),
                    found: coeff_modulus.to_string(),
                });
            }
            coeffs_cxx.push(coeff.value().to_string());
        }
        Ok(DCRTPoly::new(ffi::DCRTPolyGenFromVec(params.get_params(), &coeffs_cxx)))
    }

    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(params.get_params(), &constant.value().to_string()))
    }

    fn const_zero(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(params.get_params(), &BigUint::ZERO.to_string()))
    }

    fn const_one(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.get_params(),
            &BigUint::from(1u32).to_string(),
        ))
    }

    fn const_minus_one(params: &Self::Params) -> Self {
        DCRTPoly::new(ffi::DCRTPolyGenFromConst(
            params.get_params(),
            &(params.modulus().as_ref() - BigUint::from(1u32)).to_string(),
        ))
    }
}

impl Add for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        DCRTPoly::new(ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly))
    }
}

impl Mul for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        DCRTPoly::new(ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly))
    }
}

impl Sub for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<'a> Sub<&'a DCRTPoly> for DCRTPoly {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self + rhs.clone().neg()
    }
}

impl Neg for DCRTPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        DCRTPoly::new(self.ptr_poly.Negate())
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
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl<'a> AddAssign<&'a DCRTPoly> for DCRTPoly {
    fn add_assign(&mut self, rhs: &'a Self) {
        let res = ffi::DCRTPolyAdd(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl MulAssign for DCRTPoly {
    fn mul_assign(&mut self, rhs: Self) {
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl<'a> MulAssign<&'a DCRTPoly> for DCRTPoly {
    fn mul_assign(&mut self, rhs: &'a Self) {
        let res = ffi::DCRTPolyMul(&rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl SubAssign for DCRTPoly {
    fn sub_assign(&mut self, rhs: Self) {
        // Negate rhs and then add it to self
        let neg_rhs = rhs.neg();
        let res = ffi::DCRTPolyAdd(&neg_rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

impl<'a> SubAssign<&'a DCRTPoly> for DCRTPoly {
    fn sub_assign(&mut self, rhs: &'a Self) {
        // Clone the reference to negate it
        let neg_rhs = rhs.clone().neg();
        let res = ffi::DCRTPolyAdd(&neg_rhs.ptr_poly, &self.ptr_poly);
        self.ptr_poly = res.into();
    }
}

#[cfg(test)]
mod tests {

    use crate::poly::PolyParams;

    use super::*;

    #[test]
    fn test_dcrtpoly_arithmetic() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let q = params.modulus();

        // todo: replace value and modulus from param
        let coeffs1 = [
            FinRingElem::new(100u32, q.clone()),
            FinRingElem::new(200u32, q.clone()),
            FinRingElem::new(300u32, q.clone()),
            FinRingElem::new(400u32, q.clone()),
        ];
        let coeffs2 = [
            FinRingElem::new(500u32, q.clone()),
            FinRingElem::new(600u32, q.clone()),
            FinRingElem::new(700u32, q.clone()),
            FinRingElem::new(800u32, q.clone()),
        ];

        // 3. Create polynomials from those coefficients.
        let poly1 = DCRTPoly::from_coeffs(&params, &coeffs1).unwrap();
        let poly2 = DCRTPoly::from_coeffs(&params, &coeffs2).unwrap();

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
        let const_poly = DCRTPoly::from_const(&params, &FinRingElem::new(123, q.clone()));
        assert_eq!(
            const_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(123, q.clone()); 1]).unwrap(),
            "from_const should produce a polynomial with all coeffs = 123"
        );
        let zero_poly = DCRTPoly::const_zero(&params);
        assert_eq!(
            zero_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(0, q.clone()); 1]).unwrap(),
            "const_zero should produce a polynomial with all coeffs = 0"
        );

        let one_poly = DCRTPoly::const_one(&params);
        assert_eq!(
            one_poly,
            DCRTPoly::from_coeffs(&params, &[FinRingElem::new(1, q); 1]).unwrap(),
            "one_poly should produce a polynomial with all coeffs = 1"
        );
    }
}
