use num_bigint::BigUint;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParamsImpl},
};
use std::fmt::Debug;
use std::sync::Arc;

pub trait PolyElemParams: Clone {
    fn modulus(&self) -> BigUint;
}

pub trait PolyParams: Clone + PolyElemParams {
    fn ring_dimension(&self) -> u32;
}
