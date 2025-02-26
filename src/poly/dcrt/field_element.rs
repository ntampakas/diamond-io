use std::ops::{Add, Mul};

use num_bigint::BigUint;

use crate::poly::PElem;

#[derive(Clone, Debug, Default, Eq, PartialEq, PartialOrd, Ord)]
pub struct FieldElement {
    value: BigUint,
    modulus: BigUint,
}

impl FieldElement {
    pub fn new<T: Into<BigUint>>(value: T, modulus: T) -> Self {
        let value = value.into();
        let modulus = modulus.into();
        let reduced_value = value % modulus.clone();
        Self { value: reduced_value, modulus }
    }

    pub fn value(&self) -> BigUint {
        self.value.clone()
    }

    pub fn modulus(&self) -> BigUint {
        self.modulus.clone()
    }
}

impl PElem for FieldElement {}

// ====== Arithmetic ======

impl Add for FieldElement {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value, self.modulus)
    }
}

impl Mul for FieldElement {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.value * rhs.value, self.modulus)
    }
}
