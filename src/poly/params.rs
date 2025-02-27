use std::sync::Arc;

use num_bigint::BigUint;
use num_traits::Num;
use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParamsImpl},
};

#[derive(Clone)]
pub struct PolyParams {
    ptr_params: Arc<UniquePtr<ILDCRTParamsImpl>>,
}

impl PolyParams {
    pub fn new(n: u32, size: u32, k_res: u32) -> Self {
        let ptr_params = ffi::GenILDCRTParamsByOrderSizeBits(2 * n, size, k_res);
        Self { ptr_params: ptr_params.into() }
    }

    pub fn get_params(&self) -> &UniquePtr<ILDCRTParamsImpl> {
        &self.ptr_params
    }

    pub fn get_modulus(&self) -> BigUint {
        let modulus = &self.ptr_params.as_ref().GetModulus();
        BigUint::from_str_radix(modulus, 10).unwrap()
    }

    pub fn get_ring_dimension(&self) -> u32 {
        let ring_dimension = &self.ptr_params.as_ref().GetRingDimension();
        *ring_dimension
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correct_params_initiation() {
        let n = 16;
        let size = 4;
        let k_res = 51;
        let _ = PolyParams::new(n, size, k_res);
    }
}
