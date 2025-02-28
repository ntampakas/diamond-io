use crate::poly::{
    dcrt::{DCRTPolyMatrix, DCRTPolyParams},
    sampler::PolyTrapdoorSampler,
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{DCRTPolyTrapdoorGen, TrapdoorOutput},
};

pub struct DCRTPolyTrapdoorSampler {
    params: DCRTPolyParams,
}

pub struct DCRTPolyTrapdoor {
    _trapdoor: UniquePtr<TrapdoorOutput>,
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = DCRTPolyTrapdoor;

    fn trapdoor(&self, base: usize) -> Self::Trapdoor {
        let trapdoor = DCRTPolyTrapdoorGen(self.params.get_params(), base as i64, false);
        Self::Trapdoor { _trapdoor: trapdoor }
    }

    fn preimage(&self, _trapdoor: &Self::Trapdoor, _target: &Self::M) -> Self::M {
        todo!()
    }
}
