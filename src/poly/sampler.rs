use openfhe::ffi::{self};

use super::{
    dcrt_matrix::DCRTPolyMatrix,
    dcrt_poly::DCRTPoly,
    matrix::{get_null_matrix, PolynomialMatrix},
    params::Params,
};

pub enum DistType {
    FinRingDist,
    GaussianDist(f64),
    BitDist,
}

// pub struct FinRingDist;
// impl DistType for FinRingDist {}

// pub struct GaussianDist {
//     pub gaussian_param: f64,
// }
// impl DistType for GaussianDist {}

// pub struct BitDist;
// impl DistType for BitDist {}

pub struct MatrixUniformSampler<const ROW: usize, const COLUMNS: usize> {
    dist_type: DistType,
}

impl<const ROW: usize, const COLUMNS: usize> MatrixUniformSampler<ROW, COLUMNS> {
    pub fn new(dist_type: DistType) -> Self {
        Self { dist_type }
    }
}

impl<const ROW: usize, const COLUMNS: usize> MatrixUniformSampler<ROW, COLUMNS> {
    pub fn sample_uniform_matrix(
        &self,
        params: &Params,
    ) -> Result<DCRTPolyMatrix<DCRTPoly, ROW, COLUMNS>, anyhow::Error> {
        let mut collect = get_null_matrix::<DCRTPoly, ROW, COLUMNS>();
        #[allow(clippy::needless_range_loop)]
        for row in 0..ROW {
            for col in 0..COLUMNS {
                let sampled_poly = self.sample_poly(params)?;
                collect[row][col] = sampled_poly;
            }
        }
        let r = DCRTPolyMatrix::<DCRTPoly, ROW, COLUMNS>::from_slice(&collect);
        Ok(r)
    }

    fn sample_poly(&self, params: &Params) -> Result<DCRTPoly, anyhow::Error> {
        let sampled_poly = match self.dist_type {
            DistType::FinRingDist => ffi::DCRTPolyGenFromDug(&params.ptr_params),
            DistType::GaussianDist(sigma) => ffi::DCRTPolyGenFromDgg(&params.ptr_params, sigma),
            DistType::BitDist => ffi::DCRTPolyGenFromBug(&params.ptr_params),
        };
        Ok(DCRTPoly::new(sampled_poly))
    }
}

// pub trait PolyHashSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D: DistType>:
//     PolyDegree<T>
// {
//     type Error: std::error::Error + Send + Sync + 'static;
//     fn sample_hash<B: AsRef<[u8]>>(
//         &self,
//         tag: B,
//         rows: usize,
//         columns: usize,
//     ) -> Result<PolyMatrix<T, P, M>, Self::Error>;

//     fn set_key(&mut self, key: &[u8]);

//     fn expose_key(&self) -> &[u8];
// }

// #[allow(clippy::type_complexity)]
// pub trait PolyTrapdoorSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D:
// DistType>:     PolyDegree<T>
// {
//     type Error: std::error::Error + Send + Sync + 'static;
//     type Trapdoor;
//     fn trapdoor<R: RngCore>(
//         &self,
//         rng: &mut R,
//     ) -> Result<(PolyMatrix<T, P, M>, Self::Trapdoor), Self::Error>;

//     fn preimage<R: RngCore>(
//         &self,
//         rng: &mut R,
//         trapdoor: &Self::Trapdoor,
//         target: &PolyMatrix<T, P, M>,
//     ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt_matrix::mult;

    #[test]
    fn test_ring_dist() {
        let params = Params::new(16, 4, 51);

        // Test FinRingDist
        let sampler = MatrixUniformSampler::<20, 5>::new(DistType::FinRingDist);
        let result1 = sampler.sample_uniform_matrix(&params);
        assert!(result1.is_ok());
        let matrix1 = result1.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(&params);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let sampler = MatrixUniformSampler::<5, 12>::new(DistType::FinRingDist);
        let result3 = sampler.sample_uniform_matrix(&params);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition - TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);

        let mult_matrix: DCRTPolyMatrix<DCRTPoly, 20, 12> =
            mult::<DCRTPoly, 5, 20, 12>(&matrix1, &matrix3, params);

        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_gaussian_dist() {
        let params = Params::new(16, 4, 51);

        // Test GaussianDist
        let sampler = MatrixUniformSampler::<20, 5>::new(DistType::GaussianDist(4.57825));
        let result = sampler.sample_uniform_matrix(&params);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(&params);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let sampler = MatrixUniformSampler::<5, 12>::new(DistType::FinRingDist);
        let result3 = sampler.sample_uniform_matrix(&params);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);

        let mult_matrix: DCRTPolyMatrix<DCRTPoly, 20, 12> =
            mult::<DCRTPoly, 5, 20, 12>(&matrix1, &matrix3, params);

        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }

    #[test]
    fn test_bit_dist() {
        let params = Params::new(16, 4, 51);

        // Test BitDist
        let sampler = MatrixUniformSampler::<20, 5>::new(DistType::BitDist);
        let result = sampler.sample_uniform_matrix(&params);
        assert!(result.is_ok());
        let matrix1 = result.unwrap();
        assert_eq!(matrix1.row_size(), 20);
        assert_eq!(matrix1.col_size(), 5);

        let result2 = sampler.sample_uniform_matrix(&params);
        assert!(result2.is_ok());
        let matrix2 = result2.unwrap();

        let sampler = MatrixUniformSampler::<5, 12>::new(DistType::FinRingDist);
        let result3 = sampler.sample_uniform_matrix(&params);
        assert!(result3.is_ok());
        let matrix3 = result3.unwrap();
        assert_eq!(matrix3.row_size(), 5);
        assert_eq!(matrix3.col_size(), 12);

        // Test matrix addition - TODO: move this to DCRTPolyMatrix tests
        let added_matrix = matrix1.clone() + matrix2;
        assert_eq!(added_matrix.row_size(), 20);
        assert_eq!(added_matrix.col_size(), 5);

        let mult_matrix: DCRTPolyMatrix<DCRTPoly, 20, 12> =
            mult::<DCRTPoly, 5, 20, 12>(&matrix1, &matrix3, params);

        assert_eq!(mult_matrix.row_size(), 20);
        assert_eq!(mult_matrix.col_size(), 12);
    }
}
