use super::PolyMatrix;
use proptest::prelude::*;

#[derive(Debug, Clone, Copy)]
pub enum DistType {
    FinRingDist,
    GaussDist { sigma: f64 },
    BitDist,
}

impl Arbitrary for DistType {
    type Parameters = ();
    type Strategy = BoxedStrategy<Self>;

    fn arbitrary_with(_args: Self::Parameters) -> Self::Strategy {
        prop_oneof![
            Just(DistType::FinRingDist),
            any::<f64>().prop_map(|sigma| DistType::GaussDist { sigma }),
            Just(DistType::BitDist)
        ]
        .boxed()
    }
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

pub trait PolyUniformSampler {
    type M: PolyMatrix;

    fn sample_uniform(&self, rows: usize, columns: usize, dist: DistType) -> Self::M;
}

pub trait PolyTrapdoorSampler {
    type M: PolyMatrix;
    type Trapdoor;
    fn trapdoor(&self) -> (Self::Trapdoor, Self::M);
    fn preimage(&self, trapdoor: &Self::Trapdoor, target: &Self::M, sigma: f64) -> Self::M;
}
