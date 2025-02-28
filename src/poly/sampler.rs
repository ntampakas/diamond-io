use super::PolyMatrix;

// pub enum DistType {
//     FinRingDist,
//     GaussDist,
//     BitDist,
// }

#[derive(Debug, Clone, Copy)]
pub enum DistType {
    FinRingDist,
    GaussDist { sigma: f64 },
    BitDist,
}

// pub struct FinRingDist;
// impl DistType for FinRingDist {}
// pub struct GaussDist {}
// impl DistType for GaussDist {}

// pub struct BitDist;
// impl DistType for BitDist {}

pub trait PolyHashSampler<K: AsRef<[u8]>> {
    // type Error: std::error::Error + Send + Sync + 'static;
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
    // type Error: std::error::Error + Send + Sync + 'static;
    type M: PolyMatrix;

    fn sample_uniform(&self, rows: usize, columns: usize, dist: DistType) -> Self::M;
}

pub trait PolyTrapdoorSampler {
    // type Error: std::error::Error + Send + Sync + 'static;
    type M: PolyMatrix;
    type Trapdoor;
    fn trapdoor(&self) -> (Self::Trapdoor, Self::M);
    fn preimage(&self, trapdoor: &Self::Trapdoor, target: &Self::M, sigma: f64) -> Self::M;
}
