use std::ops::{Add, Mul};

use crate::poly::PElem;

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

impl PElem for FieldElement {}

// ====== Arithmetic ======

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
