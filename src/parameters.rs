use phantom_zone_math::{
    prelude::{ElemFrom, Modulus, Prime},
    ring::{PrimeRing, RingOps},
};

use crate::utils::empty_vector_ring;

/// Parameters for the BGG+ RLWE attribute encoding
///
/// # Fields
///
/// * `ell`: number of attributes
/// * `m`: k + 2, where k is the number of bits in the modulus
/// * `ring`: RLWE ring associated to the parameters
/// * `g`: gadget vector, which each element is a constant polynomial and there are m of them (m - 2 of them are non-zero)
#[derive(Debug, Clone)]
pub struct Parameters {
    ell: usize,
    m: usize,
    ring: PrimeRing,
    g: Vec<Vec<u64>>,
}

impl Parameters {
    /// Initialize the parameters for the BGG+ RLWE attribute encoding
    ///
    /// # Arguments
    ///
    /// * `log_ring_size`: log2 of ring size
    /// * `k`: number of bits in the ring modulus (q)
    /// * `ell`: number of attributes
    pub fn new(log_ring_size: usize, k: usize, ell: usize) -> Self {
        let q: Modulus = Prime::gen(k, log_ring_size + 1).into();
        let ring_size = 1 << log_ring_size;
        let k_ = (q.as_f64()).log2().ceil() as usize; // actual number of bits in the modulus after q is chosen
        let ring = <PrimeRing as RingOps>::new(q, ring_size);
        let m = k_ + 2;
        let g = init_gadget_vector(&ring, m);
        Self { ell, m, ring, g }
    }

    pub fn ring(&self) -> &PrimeRing {
        &self.ring
    }

    pub fn m(&self) -> &usize {
        &self.m
    }

    pub fn ell(&self) -> &usize {
        &self.ell
    }

    pub fn g(&self) -> &Vec<Vec<u64>> {
        &self.g
    }
}

/// Initialize the gadget vector `g` for the BGG+ RLWE attribute encoding
///
/// `g = [2^0, 2^1, ..., 2^(k-1), 0, 0]` where each element is a constant polynomial in the ring
pub fn init_gadget_vector(ring: &PrimeRing, m: usize) -> Vec<Vec<u64>> {
    let mut g = empty_vector_ring(ring, m);

    for i in 0..m - 2 {
        g[i][0] = ring.elem_from(2u64.pow(i as u32));
    }
    g
}
