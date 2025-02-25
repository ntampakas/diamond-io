use std::{
    fmt::Debug,
    ops::{Add, Mul},
};

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, PartialOrd, Ord)]
pub struct FieldElement {
    value: u64,
    modulus: u64,
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

impl super::PElem for FieldElement {}
