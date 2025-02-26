use std::marker::PhantomData;

use openfhe::ffi;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix},
    PolyParams,
};

use super::{MatrixUniformSamplerTrait, Polynomial, PolynomialMatrix};

pub enum DistType {
    FinRingDist,
    GaussianDist(f64),
    BitDist,
}

pub struct MatrixUniformSampler<P, M>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
{
    dist_type: DistType,
    params: P::Params,
    _phantom_m: PhantomData<M>,
}

impl MatrixUniformSamplerTrait<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
    for MatrixUniformSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
{
    type Error = std::io::Error;

    fn sample_uniform(
        &self,
        nrow: usize,
        ncol: usize,
    ) -> Result<DCRTPolyMatrix<DCRTPoly>, Self::Error> {
        let mut c: Vec<Vec<DCRTPoly>> = vec![vec![DCRTPoly::const_zero(&self.params); ncol]; nrow];
        for row in 0..nrow {
            for col in 0..ncol {
                c[row][col] = self.sample_poly(&self.params);
            }
        }
        Ok(DCRTPolyMatrix::<DCRTPoly>::from_poly_vec(&self.params, c))
    }
}

impl MatrixUniformSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>> {
    pub fn new(dist_type: DistType, params: PolyParams) -> Self {
        Self { dist_type, params, _phantom_m: PhantomData }
    }

    fn sample_poly(&self, params: &PolyParams) -> DCRTPoly {
        let sampled_poly = match self.dist_type {
            DistType::FinRingDist => ffi::DCRTPolyGenFromDug(&params.ptr_params),
            DistType::GaussianDist(sigma) => ffi::DCRTPolyGenFromDgg(&params.ptr_params, sigma),
            DistType::BitDist => ffi::DCRTPolyGenFromBug(&params.ptr_params),
        };
        DCRTPoly::new(sampled_poly)
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
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result1 = sampler.sample_uniform(20, 5);
        assert!(result1.is_ok());
        let matrix1 = result1.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform(5, 12);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition - TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_gaussian_dist() {
        let params = PolyParams::new(16, 4, 51);

        // Test GaussianDist
        let sampler = MatrixUniformSampler::new(DistType::GaussianDist(4.57825), params);
        let result = sampler.sample_uniform(20, 5);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();
        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform(5, 12);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_bit_dist() {
        let params = PolyParams::new(16, 4, 51);

        // Test BitDist
        let sampler = MatrixUniformSampler::new(DistType::BitDist, params);
        let result = sampler.sample_uniform(20, 5);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform(5, 12);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition - TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);
        let mult_matrix = matrix1 * matrix3;
        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }
}
