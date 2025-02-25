use openfhe::ffi::{self};

use super::{
    dcrt_matrix::DCRTPolyMatrix, dcrt_poly::DCRTPoly, matrix::PolynomialMatrix, params::Params,
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

pub struct MatrixUniformSampler<const ROWS: usize, const COLUMNS: usize> {
    dist_type: DistType,
}

impl<const ROWS: usize, const COLUMNS: usize> MatrixUniformSampler<ROWS, COLUMNS> {
    pub fn new(dist_type: DistType) -> Self {
        Self { dist_type }
    }
}

impl<const ROWS: usize, const COLUMNS: usize> MatrixUniformSampler<ROWS, COLUMNS> {
    pub fn sample_uniform(
        self,
        params: &Params,
    ) -> Result<DCRTPolyMatrix<DCRTPoly, ROWS, COLUMNS>, anyhow::Error> {
        let mut collect = vec![];
        for _ in 0..ROWS {
            for _ in 0..COLUMNS {
                let sampled_poly = self.sample(params)?;
                collect.push(sampled_poly);
            }
        }
        let r = DCRTPolyMatrix::<DCRTPoly, ROWS, COLUMNS>::from_slice(&collect);
        Ok(r)
    }

    fn sample(&self, params: &Params) -> Result<DCRTPoly, anyhow::Error> {
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

    #[test]
    fn test_uniform_sampler() {
        let params = Params::new(16, 4, 51);
        let sampler = MatrixUniformSampler::<1, 3>::new(DistType::FinRingDist);
        let result = sampler.sample_uniform(&params);
        assert!(result.is_ok());

        let sampler = MatrixUniformSampler::<1, 3>::new(DistType::GaussianDist(4.57825));
        let result = sampler.sample_uniform(&params);
        assert!(result.is_ok());

        let sampler = MatrixUniformSampler::<1, 3>::new(DistType::BitDist);
        let result = sampler.sample_uniform(&params);
        assert!(result.is_ok());
    }
}
