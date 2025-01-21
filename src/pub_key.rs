use crate::{parameters::Parameters, utils::empty_matrix_ring};
use phantom_zone_math::{prelude::Sampler, ring::RingOps};
use rand::thread_rng;

/// A public key in the BGG+ RLWE encoding scheme
///
/// # Fields
///
/// * `b`: The public key matrix of size `(ell + 1) x m` where each slot is a ring element
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
        let mut b = empty_matrix_ring(ring, ell + 1, m);
        for i in 0..(ell + 1) {
            for j in 0..m {
                b[i][j] = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
            }
        }
        Self { b, params }
    }

    pub fn b(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.b
    }

    pub fn params(&self) -> &Parameters {
        &self.params
    }
}
