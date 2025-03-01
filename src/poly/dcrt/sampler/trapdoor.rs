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
    ffi::{
        DCRTPolyGaussSamp, DCRTPolyTrapdoorGen, DCRTTrapdoorImpl, GetMatrixElement,
        RLWETrapdoorPair,
    },
};

pub struct DCRTPolyTrapdoor {
    // TODO: this struct should disappear and type Trapdoor should be Arc<UniquePtr<RLWETrapdoorPair>>
    trapdoor_output: Arc<UniquePtr<DCRTTrapdoorImpl>>,
    _trapdoor_pair: Arc<UniquePtr<RLWETrapdoorPair>>,
}

pub struct DCRTPolyTrapdoorSampler {
    params: DCRTPolyParams,
    base: usize,
    sigma: f64,
}

impl DCRTPolyTrapdoorSampler {
    pub fn new(params: DCRTPolyParams, base: usize, sigma: f64) -> Self {
        Self { params, base, sigma }
    }
}

impl PolyTrapdoorSampler for DCRTPolyTrapdoorSampler {
    type M = DCRTPolyMatrix;
    type Trapdoor = DCRTPolyTrapdoor;

    fn trapdoor(&self) -> (Self::Trapdoor, Self::M) {
        let trapdoor_output =
            DCRTPolyTrapdoorGen(self.params.get_params(), self.base as i64, false);
        let trapdoor_pair = trapdoor_output.GetTrapdoorPtr();
        let ncol = ceil_log2(&self.params.modulus()) + 2;

        let mut matrix_inner = Vec::with_capacity(1);
        let mut row = Vec::with_capacity(ncol);
        for i in 0..ncol {
            let poly = trapdoor_output.GetPolyAtIndex(i); // TODO: replace this functionality
            let dcrt_poly = DCRTPoly::new(poly);
            row.push(dcrt_poly);
        }
        matrix_inner.push(row);
        let row_matrix = DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner);
        let trapdoor = DCRTPolyTrapdoor {
            trapdoor_output: trapdoor_output.into(),
            _trapdoor_pair: trapdoor_pair.into(),
        };

        (trapdoor, row_matrix)
    }

    fn preimage(&self, trapdoor: &Self::Trapdoor, target: &Self::M) -> Self::M {
        let target_poly = target.entry(0, 0).clone(); // TODO: target must be a matrix
        let n = self.params.ring_dimension();
        let k = ceil_log2(&self.params.modulus());

        let preimage = DCRTPolyGaussSamp(
            n as usize,
            k,
            &trapdoor.trapdoor_output,
            target_poly.get_poly(),
            2,
            self.sigma,
        );

        let nrow = ceil_log2(&self.params.modulus()) + 2;

        let mut matrix_inner = Vec::with_capacity(nrow);
        for i in 0..nrow {
            let poly = GetMatrixElement(&preimage, i, 0); // TODO: replace this functionality
            let dcrt_poly = DCRTPoly::new(poly);
            matrix_inner.push(vec![dcrt_poly]);
        }

        DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;

    #[test]
    fn test_trapdoor_generation() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let base = 2;
        let sigma = 4.57825;
        let sampler = DCRTPolyTrapdoorSampler::new(params, base, sigma);

        let (_trapdoor, public_matrix) = sampler.trapdoor();

        let expected_cols = ceil_log2(&sampler.params.modulus()) + 2;

        // Check dimensions of the public matrix
        assert_eq!(public_matrix.row_size(), 1, "Public matrix should have 1 row");
        assert_eq!(
            public_matrix.col_size(),
            expected_cols,
            "Public matrix should have ceil_log2(q) + 2 (m) columns"
        );

        // Verify that all entries in the matrix are valid DCRTPolys
        for i in 0..public_matrix.row_size() {
            for j in 0..public_matrix.col_size() {
                let poly = public_matrix.entry(i, j);
                assert!(!poly.get_poly().is_null(), "Matrix entry should be a valid DCRTPoly");
            }
        }
    }

    #[test]
    fn test_preimage_generation() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let base = 2;
        let sigma = 4.57825;
        let sampler = DCRTPolyTrapdoorSampler::new(params.clone(), base, sigma);
        let (trapdoor, public_matrix) = sampler.trapdoor();

        // create a target matrix with 1 row and 1 column
        let target =
            DCRTPolyMatrix::from_poly_vec(&params, vec![vec![public_matrix.entry(0, 0).clone()]]);
        let preimage = sampler.preimage(&trapdoor, &target);

        // Public matrix * preimage should be equal to target
        let product = public_matrix * &preimage;
        assert_eq!(product, target, "Product of public matrix and preimage should equal target");
    }
}
