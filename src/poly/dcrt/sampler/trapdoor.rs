use std::sync::Arc;

use crate::poly::{
    dcrt::{DCRTPolyMatrix, DCRTPolyParams},
    sampler::PolyTrapdoorSampler,
};

use openfhe::{cxx::UniquePtr, ffi::RLWETrapdoorPair};

pub struct DCRTPolyTrapdoorSampler {
    _params: DCRTPolyParams,
    _base: usize,
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = Arc<UniquePtr<RLWETrapdoorPair>>;

    fn trapdoor(&self) -> (Self::Trapdoor, Self::M) {
        todo!()
        // let trapdoor = DCRTPolyTrapdoorGen(self.params.get_params(), self.base as i64, false); // TODO: check if we need balanced to be true
        // trapdoor.into()
    }

    fn preimage(&self, _trapdoor: &Self::Trapdoor, _target: &Self::M, _sigma: f64) -> Self::M {
        todo!()
        // let n_row = target.row_size();
        // let n_col = target.col_size();
        // let mut preimages = Vec::with_capacity(n_row);
        // for i in 0..n_row {
        //     let mut row_preimages = Vec::with_capacity(n_col);
        //     for j in 0..n_col {
        //         let target_poly = target.entry(i, j).clone();
        //         let preimage =
        //             DCRTPolyGaussSamp(12, 5, trapdoor.get_trapdoor(), &target_poly.get_poly(), 10);
        //         row_preimages.push(preimage);
        //     }
        //     preimages.push(row_preimages);
        // }
        // Self::M::from_poly_vec(&self.params, preimages)
    }
}
