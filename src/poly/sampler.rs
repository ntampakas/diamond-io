use super::PolyMatrix;

#[derive(Debug, Clone, Copy)]
pub enum DistType {
    FinRingDist,
    GaussDist { sigma: f64 },
    BitDist,
}

pub trait PolyHashSampler<K: AsRef<[u8]>> {
    type M: PolyMatrix;

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        rows: usize,
        columns: usize,
        dist: DistType,
    ) -> Self::M;

    fn set_key(&mut self, key: K);

    fn expose_key(&self) -> &[u8];
}

pub trait PolyUniformSample {
    type M: PolyMatrix;

    fn sample_uniform(&self, rows: usize, columns: usize, dist: DistType) -> Self::M;
}

pub trait PolyTrapdoorSampler {
    type M: PolyMatrix;
    type Trapdoor;
    fn trapdoor(&self) -> (Self::M, Self::Trapdoor);

    fn preimage(&self, trapdoor: &Self::Trapdoor, target: &Self::M) -> Self::M;
}
