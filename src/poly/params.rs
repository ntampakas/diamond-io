use std::{fmt::Debug, sync::Arc};

use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParams},
};

#[derive(Clone)]
pub struct Params {
    pub ptr_params: Arc<UniquePtr<ILDCRTParams>>,
    n: u32,
    size: u32,
    k_res: u32,
}

impl Debug for Params {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Params")
            .field("n", &self.n) // Added field n
            .field("size", &self.size) // Added field size
            .field("k_res", &self.k_res) // Added field k_res
            .finish()
    }
}
impl Params {
    pub fn new(n: u32, size: u32, k_res: u32) -> Self {
        let ptr_params = ffi::GenILDCRTParamsByOrderSizeBits(2 * n, size, k_res);
        Self { ptr_params: ptr_params.into(), n, size, k_res }
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
        let _ = Params::new(n, size, k_res);
    }
}
