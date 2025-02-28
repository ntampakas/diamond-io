use crate::poly::params::PolyParams;
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

pub trait PolyElem:
    Sized
    + Debug
    + Eq
    + Ord
    + Send
    + Sync
    + Clone
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
{
    type Params;
    fn zero(params: &Self::Params) -> Self;
    fn one(params: &Self::Params) -> Self;
    fn minus_one(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self) -> bool;
}

/// Describes the common interface polynomials
pub trait Poly:
    Sized
    + Clone
    + Debug
    + PartialEq
    + Eq
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + AddAssign
    + SubAssign
    + MulAssign
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
{
    type Elem: PolyElem;
    type Params: PolyParams;
    fn from_coeffs(params: &Self::Params, coeffs: &[Self::Elem]) -> Self;
    fn from_const(params: &Self::Params, constant: &Self::Elem) -> Self;
    fn coeffs(&self) -> Vec<Self::Elem>;
    fn const_zero(params: &Self::Params) -> Self;
    fn const_one(params: &Self::Params) -> Self;
    fn const_minus_one(params: &Self::Params) -> Self;
    fn extract_highest_bits(&self) -> Vec<bool> {
        let mut bits = Vec::new();
        for elem in self.coeffs() {
            bits.push(elem.extract_highest_bits());
        }
        bits
    }
}
