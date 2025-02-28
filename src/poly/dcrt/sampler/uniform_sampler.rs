use openfhe::ffi;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
    sampler::{DistType, PolyUniformSample},
    Poly, PolyMatrix,
};

pub struct DCRTPolyUniformSampler {
    params: DCRTPolyParams,
}

impl DCRTPolyUniformSampler {
    pub fn new(params: DCRTPolyParams) -> Self {
        Self { params }
    }

    fn sample_poly(&self, params: &DCRTPolyParams, dist: &DistType) -> DCRTPoly {
        let sampled_poly = match dist {
            DistType::FinRingDist => ffi::DCRTPolyGenFromDug(params.get_params()),
            DistType::GaussDist { sigma } => ffi::DCRTPolyGenFromDgg(params.get_params(), *sigma),
            DistType::BitDist => ffi::DCRTPolyGenFromBug(params.get_params()),
        };
        DCRTPoly::new(sampled_poly)
    }
}

impl PolyUniformSample for DCRTPolyUniformSampler {
    type M = DCRTPolyMatrix;

    fn sample_uniform(&self, rows: usize, columns: usize, dist: DistType) -> Self::M {
        let mut c: Vec<Vec<DCRTPoly>> =
            vec![vec![DCRTPoly::const_zero(&self.params); columns]; rows];
        for row in 0..rows {
            for col in 0..columns {
                let sampled_poly = self.sample_poly(&self.params, &dist);
                if sampled_poly.get_poly().is_null() {
                    panic!("Attempted to dereference a null pointer");
                }
                c[row][col] = sampled_poly;
            }
        }
        DCRTPolyMatrix::from_poly_vec(&self.params, c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ring_dist() {
        let params = DCRTPolyParams::new(16, 4, 51);

        // Test FinRingDist
        let sampler = DCRTPolyUniformSampler::new(params.clone());
        let matrix1 = sampler.sample_uniform(20, 5, DistType::FinRingDist);
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let matrix2 = sampler.sample_uniform(20, 5, DistType::FinRingDist);

        let sampler2 = DCRTPolyUniformSampler::new(params);
        let matrix3 = sampler2.sample_uniform(5, 12, DistType::FinRingDist);
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

    #[test]
    fn test_gaussian_dist() {
        let params = DCRTPolyParams::new(16, 4, 51);

        // Test GaussianDist
        let sampler = DCRTPolyUniformSampler::new(params.clone());
        let matrix1 = sampler.sample_uniform(20, 5, DistType::GaussDist { sigma: 4.57825 });
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let matrix2 = sampler.sample_uniform(20, 5, DistType::GaussDist { sigma: 4.57825 });

        let sampler2 = DCRTPolyUniformSampler::new(params);
        let matrix3 = sampler2.sample_uniform(5, 12, DistType::FinRingDist);
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

    #[test]
    fn test_bit_dist() {
        let params = DCRTPolyParams::new(16, 4, 51);

        // Test BitDist
        let sampler = DCRTPolyUniformSampler::new(params.clone());
        let matrix1 = sampler.sample_uniform(20, 5, DistType::BitDist);
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);
        // [TODO] Test the norm of each coefficient of polynomials in the matrix.

        let matrix2 = sampler.sample_uniform(20, 5, DistType::BitDist);

        let sampler2 = DCRTPolyUniformSampler::new(params);
        let matrix3 = sampler2.sample_uniform(5, 12, DistType::FinRingDist);
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
