use std::{
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::Arc,
};

use crate::poly::PolyElem;
use num_bigint::{BigInt, BigUint};

#[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
pub struct FinRing {
    value: BigUint,
    modulus: Arc<BigUint>,
}

impl FinRing {
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

impl PolyElem for FinRing {
    type Modulus = Arc<BigUint>;
    fn zero(modulus: &Self::Modulus) -> Self {
        Self::new(0, modulus.clone())
    }

    fn one(modulus: &Self::Modulus) -> Self {
        Self::new(1, modulus.clone())
    }

    fn minus_one(modulus: &Self::Modulus) -> Self {
        let max_minus_one = modulus.as_ref() - &BigUint::from(1u8);
        Self::new(max_minus_one, modulus.clone())
    }

    fn extract_highest_bits(&self) -> bool {
        self.value < self.modulus() / 2u8
    }
}

// ====== Arithmetic ======

impl Add for FinRing {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        self + &rhs
    }
}

impl<'a> Add<&'a FinRing> for FinRing {
    type Output = Self;

    fn add(self, rhs: &'a FinRing) -> Self::Output {
        Self::new(self.value + &rhs.value, self.modulus)
    }
}

impl AddAssign for FinRing {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value + rhs.value, self.modulus.clone());
    }
}

impl<'a> AddAssign<&'a FinRing> for FinRing {
    fn add_assign(&mut self, rhs: &'a FinRing) {
        *self = Self::new(&self.value + &rhs.value, self.modulus.clone());
    }
}

impl Mul for FinRing {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        self * &rhs
    }
}

impl<'a> Mul<&'a FinRing> for FinRing {
    type Output = Self;

    fn mul(self, rhs: &'a FinRing) -> Self::Output {
        Self::new(self.value * &rhs.value, self.modulus)
    }
}

impl MulAssign for FinRing {
    fn mul_assign(&mut self, rhs: Self) {
        *self = Self::new(&self.value * rhs.value, self.modulus.clone());
    }
}

impl<'a> MulAssign<&'a FinRing> for FinRing {
    fn mul_assign(&mut self, rhs: &'a FinRing) {
        *self = Self::new(&self.value * &rhs.value, self.modulus.clone());
    }
}

impl Sub for FinRing {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl<'a> Sub<&'a FinRing> for FinRing {
    type Output = Self;

    fn sub(self, rhs: &'a FinRing) -> Self::Output {
        self + (&-rhs.clone())
    }
}

impl SubAssign for FinRing {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs;
    }
}

impl<'a> SubAssign<&'a FinRing> for FinRing {
    fn sub_assign(&mut self, rhs: &'a FinRing) {
        *self = self.clone() - rhs;
    }
}

impl Neg for FinRing {
    type Output = Self;

    fn neg(self) -> Self::Output {
        if self.value == BigUint::from(0u8) {
            return self;
        }
        Self::new(self.modulus.as_ref() - &self.value, self.modulus)
    }
}
