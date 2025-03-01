use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

use crate::poly::PolyElem;
use num_bigint::{BigInt, BigUint};

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRingElem {
    value: BigUint,
    modulus: Arc<BigUint>,
}

impl FinRingElem {
    pub fn new<V: Into<BigInt>>(value: V, modulus: Arc<BigUint>) -> Self {
        let value = value.into().to_biguint().unwrap();
        let reduced_value =
            if &value < modulus.as_ref() { value.clone() } else { value % modulus.as_ref() };
        Self { value: reduced_value, modulus }
    }

    pub fn value(&self) -> &BigUint {
        &self.value
    }

    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }
}

impl PolyElem for FinRingElem {
    type Params = Arc<BigUint>;

    fn zero(params: &Self::Params) -> Self {
        Self::new(0, params.clone())
    }

    fn one(params: &Self::Params) -> Self {
        Self::new(1, params.clone())
    }

    fn minus_one(params: &Self::Params) -> Self {
        let max_minus_one = params.as_ref() - &BigUint::from(1u8);
        Self::new(max_minus_one, params.clone())
    }

    fn extract_highest_bits(&self) -> bool {
        self.value < self.modulus() / 2u8
    }
}

// ====== Arithmetic ======

impl Add for FinRingElem {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn add(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.value + &rhs.value, self.modulus)
    }
}

impl AddAssign for FinRingElem {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value + rhs.value, self.modulus.clone());
    }
}

impl<'a> AddAssign<&'a FinRingElem> for FinRingElem {
    fn add_assign(&mut self, rhs: &'a Self) {
        *self = Self::new(&self.value + &rhs.value, self.modulus.clone());
    }
}

impl Mul for FinRingElem {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn mul(self, rhs: &'a Self) -> Self::Output {
        Self::new(self.value * &rhs.value, self.modulus)
    }
}

impl MulAssign for FinRingElem {
    fn mul_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value * rhs.value, self.modulus.clone());
    }
}

impl<'a> MulAssign<&'a FinRingElem> for FinRingElem {
    fn mul_assign(&mut self, rhs: &'a Self) {
        *self = Self::new(&self.value * &rhs.value, self.modulus.clone());
    }
}

impl Sub for FinRingElem {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<'a> Sub<&'a FinRingElem> for FinRingElem {
    type Output = Self;

    fn sub(self, rhs: &'a Self) -> Self::Output {
        self + (&-rhs.clone())
    }
}

impl SubAssign for FinRingElem {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl<'a> SubAssign<&'a FinRingElem> for FinRingElem {
    fn sub_assign(&mut self, rhs: &'a Self) {
        *self = self.clone() - rhs;
    }
}

impl Neg for FinRingElem {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.value == BigUint::from(0u8) {
            return self;
        }
        Self::new(self.modulus.as_ref() - &self.value, self.modulus)
    }
}
