use crate::operations::{bit_decompose, poly_add};
use crate::Parameters;
use phantom_zone_math::{
    prelude::{ModulusOps, Sampler},
    ring::RingOps,
};
use rand::thread_rng;

/// A public key in the BGG+ RLWE encoding scheme
///
/// # Fields
///
/// * `b`: The public key matrix of size `(ell + 1) x m` where each slot is a polynomial
/// * `params`: The system parameters
#[derive(Debug, Clone)]
pub struct PublicKey {
    b: Vec<Vec<Vec<u64>>>,
    params: Parameters,
}

impl PublicKey {
    pub fn new(params: Parameters) -> Self {
        let mut rng = thread_rng();
        let ring = params.ring();
        let m = *params.m();
        let ell = *params.ell();
        let mut b = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; ell + 1];
        for i in 0..(ell + 1) {
            for j in 0..m {
                b[i][j] = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
            }
        }
        Self { b, params }
    }

    /// Perform a gate addition over the public key components at indices `idx_a` and `idx_b` and return the result
    pub fn add_gate(&self, idx_a: usize, idx_b: usize) -> Vec<Vec<u64>> {
        let ring = self.params.ring();
        let m = *self.params.m();
        let mut out = vec![vec![ring.zero(); ring.ring_size()]; m];
        for i in 0..m {
            out[i] = poly_add(&ring, &self.b[idx_a][i], &self.b[idx_b][i]);
        }
        out
    }

    /// Perform a gate multiplication over the public key components at indices `idx_a` and `idx_b` and return the result
    pub fn mul_gate(&self, idx_a: usize, idx_b: usize) -> Vec<Vec<u64>> {
        let ring = self.params.ring();
        let m = *self.params.m();
        let mut out = vec![vec![ring.zero(); ring.ring_size()]; m];

        // Compute minus_b_idx_a by multiplying each coefficient by -1
        let mut minus_b_idx_a = vec![vec![ring.zero(); ring.ring_size()]; m];
        for i in 0..m {
            for j in 0..ring.ring_size() {
                // To get -1 * coefficient in the ring, we subtract the coefficient from 0
                minus_b_idx_a[i][j] = ring.sub(&ring.zero(), &self.b[idx_a][i][j]);
            }
        }

        let tau = bit_decompose(self.params(), &minus_b_idx_a);

        // Compute out = b_idx_b * TAU
        for i in 0..m {
            for h in 0..m {
                let mut scratch = ring.allocate_scratch(1, 2, 0);
                let mut scratch = scratch.borrow_mut();
                let product = ring.take_poly(&mut scratch);

                // Multiply b_idx_b[h] by tau[h][i]
                ring.poly_mul(product, &self.b[idx_b][h], &tau[h][i], scratch.reborrow());

                out[i] = poly_add(ring, &out[i], &product.to_vec());
            }
        }
        out
    }

    pub fn b(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.b
    }

    pub fn params(&self) -> &Parameters {
        &self.params
    }
}
