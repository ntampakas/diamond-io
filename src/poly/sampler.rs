use rand::RngCore;

use super::{
    matrix::{PolyMatrix, PolyMatrixOps},
    PElem, Poly, PolyElemOps, PolyGaussOps, PolyOps,
};

pub trait PolyUniformSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    type Error: std::error::Error;
    fn sample<R: RngCore>(
        &self,
        rng: &mut R,
        rows: usize,
        columns: usize,
    ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
}

pub trait PolyHashSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    type Error: std::error::Error;
    fn sample<B: AsRef<[u8]>>(
        &self,
        key: B,
        tag: B,
        rows: usize,
        columns: usize,
    ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
}

pub trait PolyTrapdoorSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>> {
    type Error: std::error::Error;
    type Trapdoor;
    fn trapdoor<R: RngCore>(
        &self,
        rng: &mut R,
    ) -> Result<(PolyMatrix<T, P, M>, Self::Trapdoor), Self::Error>;

    fn preimage<
        R: RngCore,
        U: PolyElemOps,
        PU: PolyOps<U>,
        MU: PolyMatrixOps<U, PU>,
        V: PolyGaussOps,
        PV: PolyOps<V>,
        MV: PolyMatrixOps<V, PV>,
    >(
        &self,
        rng: &mut R,
        trapdoor: &Self::Trapdoor,
        target: &PolyMatrix<U, PU, MU>,
    ) -> Result<PolyMatrix<V, PV, MV>, Self::Error>;
}
