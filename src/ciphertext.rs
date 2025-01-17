use crate::{operations::bit_decompose, poly_add, poly_sub, Parameters, PublicKey};
use phantom_zone_crypto::util::distribution::NoiseDistribution;
use phantom_zone_math::{
    prelude::{ElemFrom, Gaussian, ModulusOps, Sampler},
    ring::RingOps,
};
use rand::{thread_rng, Rng};

/// A ciphertext in the BGG+ RLWE encoding scheme
///
/// # Fields
///
/// * `inner`: The inner component matrix of size (ell + 1) × m where the rows are equal to the vectors G+b_0, x_1G+b_1, ..., x_ellG+b_ell
/// * `secret`: The secret vector which is a polynomial in the ring
/// * `error`: The error matrix of size (ell + 1) × m where each row is equal to e_aS_0, e_aS_1, ..., e_aS_ell
/// * `params`: The system parameters
/// * `x`: The attribute vector
/// * `pub_key`: The public key
#[derive(Clone, Debug)]
pub struct Ciphertext {
    inner: Vec<Vec<Vec<u64>>>,
    secret: Vec<u64>,
    error: Vec<Vec<Vec<u64>>>,
    params: Parameters,
    x: Vec<u64>,
}

impl Ciphertext {
    /// Create a new ciphertext by encoding an attribute vector
    ///
    /// # Arguments
    /// * `public_key`: The public key matrix
    /// * `params`: System parameters
    /// * `x`: Attribute vector to encode
    pub fn new(public_key: &PublicKey, params: &Parameters, x: &Vec<u64>) -> Self {
        assert!(x.len() == params.ell + 1);

        let mut rng = thread_rng();
        let ring = &params.ring;
        let s = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
        let mut ct_inner = public_key.b().clone();

        // Generate error vectors
        let mut err_a = vec![vec![ring.zero(); ring.ring_size()]; params.m];
        let gaussian: NoiseDistribution = Gaussian(3.19).into();

        for i in 0..params.m {
            ring.sample_into::<i64>(&mut err_a[i], gaussian, rng.clone());
        }

        // Initialize error matrix
        let mut error = vec![vec![vec![ring.zero(); ring.ring_size()]; params.m]; params.ell + 1];

        for i in 0..params.ell + 1 {
            for si in 0..params.m {
                for sj in 0..params.m {
                    let random_bit = if rng.gen_bool(0.5) { 1 } else { -1 };
                    if random_bit == 1 {
                        error[i][si] = poly_add(ring, &error[i][si], &err_a[sj]);
                    }
                    if random_bit == -1 {
                        error[i][si] = poly_sub(ring, &error[i][si], &err_a[sj]);
                    }
                }
            }

            // Add gadget vector to ct_inner if x[i] = 1
            for j in 0..params.m {
                if x[i] == 1 {
                    ct_inner[i][j] = poly_add(ring, &ct_inner[i][j], &params.g[j]);
                }
            }
        }

        Self {
            inner: ct_inner,
            secret: s,
            error,
            params: params.clone(),
            x: x.clone(),
        }
    }

    /// Compute the ct_full as inner * secret + error
    pub fn compute_ct_full(&self) -> Vec<Vec<Vec<u64>>> {
        let ring = &self.params.ring;
        let mut ct_full = self.inner.clone();

        for i in 0..self.inner.len() {
            for j in 0..self.inner[0].len() {
                let mut scratch = ring.allocate_scratch(1, 2, 0);
                let mut scratch = scratch.borrow_mut();
                let c = ring.take_poly(&mut scratch);

                // Compute inner * secret
                ring.poly_mul(c, &self.inner[i][j], &self.secret, scratch.reborrow());

                // Add error
                ct_full[i][j] = super::poly_add(ring, &c.to_vec(), &self.error[i][j]);
            }
        }

        ct_full
    }

    /// Perform an add gate over the ciphertext components at indices `idx_a` and `idx_b` and return matrix `out`
    /// that satisfies the homorphism property described in [DDP+17] at page 22
    pub fn add_gate(&self, idx_a: usize, idx_b: usize) -> Vec<Vec<Vec<u64>>> {
        let ring = &self.params.ring;
        let m = self.params.m;
        let mut out = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; 2 * m];

        // Fill both identity matrices at once
        for i in 0..m {
            // Set constant polynomial 1 on diagonals of both identity matrices
            out[i][i][0] = ring.elem_from(1u64);
            out[i + m][i][0] = ring.elem_from(1u64);
        }

        out
    }

    /// Perform a mul gate over the ciphertext components at indices `idx_a` and `idx_b` and return matrix `out`
    /// that satisfies the homorphism property described in [DDP+17] at page 22
    pub fn mul_gate(&self, pub_key: &PublicKey, idx_a: usize, idx_b: usize) -> Vec<Vec<Vec<u64>>> {
        let ring = &self.params.ring;
        let m = self.params.m;

        let mut out = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; 2 * m];

        // First matrix: Identity matrix scaled by x_idx_b
        let mut x = vec![ring.zero(); ring.ring_size()];
        let x_idx_b = self.x[idx_b];
        x[0] = ring.elem_from(x_idx_b);
        for i in 0..m {
            out[i][i] = x.clone();
        }

        // Second matrix: Tau(-b_idx_a)
        // First compute -b_idx_a
        let b_idx_a = &pub_key.b()[idx_a];
        let mut minus_b_idx_a = vec![vec![ring.zero(); ring.ring_size()]; m];
        for i in 0..m {
            for j in 0..ring.ring_size() {
                minus_b_idx_a[i][j] = ring.sub(&ring.zero(), &b_idx_a[i][j]);
            }
        }

        // Compute tau of -b_idx_a
        let tau = bit_decompose(&self.params, &minus_b_idx_a);

        // Copy tau into the bottom half of out
        for h in 0..m {
            for i in 0..m {
                out[h + m][i] = tau[h][i].clone();
            }
        }

        out
    }

    pub fn inner(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.inner
    }
}
