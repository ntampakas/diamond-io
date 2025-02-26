use std::sync::Arc;

use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParamsImpl},
};

#[derive(Clone)]
pub struct PolyParams {
    pub ptr_params: Arc<UniquePtr<ILDCRTParamsImpl>>, //TODO: add getter for params
}

impl PolyParams {
    pub fn new(n: u32, size: u32, k_res: u32) -> Self {
        let ptr_params = ffi::GenILDCRTParamsByOrderSizeBits(2 * n, size, k_res);
        Self { ptr_params: ptr_params.into() }
    }

    pub fn get_modulus(&self) -> u64 {
        let modulus = &self.ptr_params.as_ref().GetModulus();
        *modulus
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
