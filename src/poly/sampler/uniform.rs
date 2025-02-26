use openfhe::ffi::{self};

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix},
    PolyParams,
};

use super::{Polynomial, PolynomialMatrix};

pub enum DistType {
    FinRingDist,
    GaussianDist(f64),
    BitDist,
}

pub struct MatrixUniformSampler {
    dist_type: DistType,
    params: PolyParams,
}

impl MatrixUniformSampler {
    pub fn new(dist_type: DistType, params: PolyParams) -> Self {
        Self { dist_type, params }
    }
}

impl MatrixUniformSampler {
    pub fn sample_uniform_matrix(
        &self,
        nrow: usize,
        ncol: usize,
    ) -> Result<DCRTPolyMatrix<DCRTPoly>, anyhow::Error> {
        let mut c: Vec<Vec<DCRTPoly>> = vec![vec![DCRTPoly::null(); ncol]; nrow];
        for row in 0..nrow {
            for col in 0..ncol {
                let sampled_poly = self.sample_poly(&self.params)?;
                c[row][col] = sampled_poly;
            }
        }
        let r = DCRTPolyMatrix::<DCRTPoly>::from_poly_vec(&self.params, c);
        Ok(r)
    }

    fn sample_poly(&self, params: &PolyParams) -> Result<DCRTPoly, anyhow::Error> {
        let sampled_poly = match self.dist_type {
            DistType::FinRingDist => ffi::DCRTPolyGenFromDug(&params.ptr_params),
            DistType::GaussianDist(sigma) => ffi::DCRTPolyGenFromDgg(&params.ptr_params, sigma),
            DistType::BitDist => ffi::DCRTPolyGenFromBug(&params.ptr_params),
        };
        Ok(DCRTPoly::new(sampled_poly))
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
        let result1 = sampler.sample_uniform_matrix(20, 5);
        assert!(result1.is_ok());
        let matrix1 = result1.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform_matrix(5, 12);
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
        let result = sampler.sample_uniform_matrix(20, 5);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();
        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform_matrix(5, 12);
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
        let result = sampler.sample_uniform_matrix(20, 5);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(20, 5);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let params = PolyParams::new(16, 4, 51);
        let sampler = MatrixUniformSampler::new(DistType::FinRingDist, params);
        let result3 = sampler.sample_uniform_matrix(5, 12);
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
