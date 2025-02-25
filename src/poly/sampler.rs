use rand::RngCore;

use super::{
    matrix::{PolyMatrix, PolyMatrixOps},
    PolyDegree, PolyElemOps, PolyOps,
};

pub trait DistType {}

pub struct FinRingDist;
impl DistType for FinRingDist {}

pub struct GaussianDist {
    pub gaussian_param: f64,
}
impl DistType for GaussianDist {}

pub struct BitDist;
impl DistType for BitDist {}

pub trait PolyUniformSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D: DistType>:
    PolyDegree<T>
{
    type Error: std::error::Error + Send + Sync + 'static;
    fn sample_uniform<R: RngCore>(
        &self,
        rng: &mut R,
        rows: usize,
        columns: usize,
    ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
}

pub trait PolyHashSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D: DistType>:
    PolyDegree<T>
{
    type Error: std::error::Error + Send + Sync + 'static;
    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        rows: usize,
        columns: usize,
    ) -> Result<PolyMatrix<T, P, M>, Self::Error>;

    fn set_key(&mut self, key: &[u8]);

    fn expose_key(&self) -> &[u8];
}

#[allow(clippy::type_complexity)]
pub trait PolyTrapdoorSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D: DistType>:
    PolyDegree<T>
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Trapdoor;
    fn trapdoor<R: RngCore>(
        &self,
        rng: &mut R,
    ) -> Result<(PolyMatrix<T, P, M>, Self::Trapdoor), Self::Error>;

    fn preimage<R: RngCore>(
        &self,
        rng: &mut R,
        trapdoor: &Self::Trapdoor,
        target: &PolyMatrix<T, P, M>,
    ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
}
