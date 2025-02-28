use std::{
    ops::{Add, Mul},
    sync::Arc,
};

use crate::poly::PolyElem;
use num_bigint::{BigInt, BigUint};

// #[derive(Clone, Debug, Eq, PartialEq, PartialOrd, Ord)]
// pub struct FinRingParams {
//     modulus: Arc<BigUint>,
// }

// impl PolyElemParams for FinRingParams {
//     fn modulus(&self) -> BigUint {
//         self.modulus.as_ref().clone()
//     }
// }

// impl FinRingParams {
//     pub fn new(modulus: Arc<BigUint>) -> Self {
//         Self { modulus }
//     }
// }

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

impl Add for FinRing {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.value + rhs.value, self.modulus)
    }
}

impl Mul for FinRing {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(self.value * rhs.value, self.modulus)
    }
}
