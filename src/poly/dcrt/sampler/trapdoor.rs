use std::sync::Arc;

use crate::{
    poly::{
        dcrt::{DCRTPoly, DCRTPolyMatrix, DCRTPolyParams},
        sampler::PolyTrapdoorSampler,
        PolyMatrix, PolyParams,
    },
    utils::ceil_log2,
};

use openfhe::{
    cxx::UniquePtr,
    ffi::{DCRTPolyTrapdoorGen, RLWETrapdoorPair},
};

pub struct DCRTPolyTrapdoorSampler {
    params: DCRTPolyParams,
    base: usize,
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = Arc<UniquePtr<RLWETrapdoorPair>>;

    fn trapdoor(&self) -> (Self::Trapdoor, Self::M) {
        let trapdoor_output =
            DCRTPolyTrapdoorGen(self.params.get_params(), self.base as i64, false);
        let trapdoor = trapdoor_output.GetTrapdoorPtr();
        let ncol = ceil_log2(&self.params.modulus()) + 2;

        let mut matrix_inner = Vec::with_capacity(1);
        let mut row = Vec::with_capacity(ncol);
        for i in 0..ncol {
            let poly = trapdoor_output.GetPolyAtIndex(i);
            let dcrt_poly = DCRTPoly::new(poly);
            row.push(dcrt_poly);
        }
        matrix_inner.push(row);
        let matrix = DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner);
        (trapdoor.into(), matrix)
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
