use num_bigint::BigUint;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParamsImpl},
};
use std::fmt::Debug;
use std::sync::Arc;

pub trait PolyElemParams: Debug + Clone {
    fn get_modulus(&self) -> BigUint;
}

pub trait PolyParams: Debug + Clone + PolyElemParams {
    fn get_ring_dimension(&self) -> u32;
}
