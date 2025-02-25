use openfhe::{
    cxx::UniquePtr,
    ffi::{self, ILDCRTParams},
};

pub struct Params {
    pub ptr_params: UniquePtr<ILDCRTParams>,
}

impl Params {
    pub fn new(n: u32, size: u32, k_res: u32) -> Self {
        let ptr_params = ffi::GenILDCRTParamsByOrderSizeBits(2 * n, size, k_res);
        Self { ptr_params }
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
