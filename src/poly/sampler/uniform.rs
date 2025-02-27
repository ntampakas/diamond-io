use std::marker::PhantomData;

use openfhe::ffi;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix},
    PolyParams,
};
use thiserror::Error;

use super::{Poly, PolynomialMatrix, UniformSampler};

#[derive(Error, Debug)]
pub enum MatrixUniformSamplerError {
    #[error("Attempted to dereference a null pointer")]
    NullPointer,
}

pub enum UniformDistType {
    FinRingDist,
    GaussianDist(f64),
    BitDist,
}

pub struct UniformSamplerImpl<P, M>
where
    P: Poly,
    M: PolynomialMatrix<P>,
{
    dist_type: UniformDistType,
    params: P::Params,
    _phantom_m: PhantomData<M>,
}

impl UniformSamplerImpl<DCRTPoly, DCRTPolyMatrix<DCRTPoly>> {
    pub fn new(dist_type: UniformDistType, params: PolyParams) -> Self {
        Self { dist_type, params, _phantom_m: PhantomData }
    }

    fn sample_poly(&self, params: &PolyParams) -> DCRTPoly {
        let sampled_poly = match self.dist_type {
            UniformDistType::FinRingDist => ffi::DCRTPolyGenFromDug(params.get_params()),
            UniformDistType::GaussianDist(sigma) => {
                ffi::DCRTPolyGenFromDgg(params.get_params(), sigma)
            }
            UniformDistType::BitDist => ffi::DCRTPolyGenFromBug(params.get_params()),
        };
        DCRTPoly::new(sampled_poly)
    }
}

impl UniformSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
    for UniformSamplerImpl<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
{
    type Error = MatrixUniformSamplerError;

    fn sample_uniform(
        &self,
        nrow: usize,
        ncol: usize,
    ) -> Result<DCRTPolyMatrix<DCRTPoly>, Self::Error> {
        let mut c: Vec<Vec<DCRTPoly>> = vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];
        for row in 0..nrow {
            for col in 0..ncol {
                let sampled_poly = self.sample_poly(&self.params);
                if sampled_poly.get_poly().is_null() {
                    return Err(MatrixUniformSamplerError::NullPointer);
                }
                c[row][col] = self.sample_poly(&self.params);
            }
        }
        Ok(DCRTPolyMatrix::<DCRTPoly>::from_poly_vec(&self.params, c))
    }
}

#[cfg(test)]
mod tests {
    use crate::poly::PolyParams;

    use super::*;

    #[test]
    fn test_ring_dist() {
        let params = PolyParams::new(16, 4, 51);

        // Test FinRingDist
        let sampler = UniformSamplerImpl::new(UniformDistType::FinRingDist, params);
        let result1 = sampler.sample_uniform(20, 5);
        assert!(result1.is_ok());
        let matrix1 = result1.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = UniformSamplerImpl::new(UniformDistType::FinRingDist, params);
        let result3 = sampler.sample_uniform(5, 12);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix arithmetic
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    // #[test] // TODO: fix this test as sometimes it throws "SIGSEGV: invalid memory reference"
    // fn test_gaussian_dist() {
    //     let params = PolyParams::new(16, 4, 51);

    //     // Test GaussianDist
    //     let sampler = MatrixUniformSampler::new(UniformDistType::GaussianDist(4.57825), params);
    //     let result = sampler.sample_uniform(20, 5);
    //     assert!(result.is_ok());
    //     let matrix1 = result.unwrap();
    //     assert_eq!(matrix1.row_size(), 20);
    //     assert_eq!(matrix1.col_size(), 5);

    //     let result2 = sampler.sample_uniform(20, 5);
    //     assert!(result2.is_ok());
    //     let matrix2 = result2.unwrap();
    //     let params = PolyParams::new(16, 4, 51);
    //     let sampler = MatrixUniformSampler::new(UniformDistType::FinRingDist, params);
    //     let result3 = sampler.sample_uniform(5, 12);
    //     assert!(result3.is_ok());
    //     let matrix3 = result3.unwrap();
    //     assert_eq!(matrix3.row_size(), 5);
    //     assert_eq!(matrix3.col_size(), 12);

    //     // Test matrix arithmetic
    //     let added_matrix = matrix1.clone() + matrix2;
    //     assert_eq!(added_matrix.row_size(), 20);
    //     assert_eq!(added_matrix.col_size(), 5);
    //     let mult_matrix = matrix1 * matrix3;
    //     assert_eq!(mult_matrix.row_size(), 20);
    //     assert_eq!(mult_matrix.col_size(), 12);
    // }

    #[test]
    fn test_bit_dist() {
        let params = PolyParams::new(16, 4, 51);

        // Test BitDist
        let sampler = UniformSamplerImpl::new(UniformDistType::BitDist, params);
        let result = sampler.sample_uniform(20, 5);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = UniformSamplerImpl::new(UniformDistType::FinRingDist, params);
        let result3 = sampler.sample_uniform(5, 12);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix arithmetic
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }
}
