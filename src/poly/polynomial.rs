use phantom_zone_math::util::serde::Serde;

use crate::poly::params::{PolyElemParams, PolyParams};
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg},
};

pub trait PolyElem: 'static + Debug + Default + Eq + Ord + Send + Sync + Serde + Add + Mul {
    type Params: PolyElemParams;

    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn minus_one(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self) -> bool;
}

/// Describes the common interface polynomials
pub trait Poly:
    Sized + Clone + Debug + PartialEq + Eq + Add + AddAssign + Mul + MulAssign + Neg
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Elem: PolyElem;
    type Params: PolyParams;
    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Result<Self, Self::Error>;
    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self;
    fn const_zero(params: &Self::Params) -> Self;
    fn const_one(params: &Self::Params) -> Self;
    fn const_minus_one(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self, poly: &Self::Poly) -> Vec<bool> {
        let mut bits = Vec::new();
        for elem in poly.coeffs() {
            bits.push(elem.extract_highest_bits());
        }
        bits
    }
}
