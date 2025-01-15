use phantom_zone_crypto::util::distribution::NoiseDistribution;
use phantom_zone_math::{
    prelude::{ElemFrom, Gaussian, ModulusOps, Sampler},
    ring::{PrimeRing, RingOps},
};
use rand::{thread_rng, Rng};

use crate::{BggRlwe, Parameters};

pub fn poly_add(ring: &PrimeRing, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    let mut c = vec![ring.zero(); a.len()];
    for i in 0..a.len() {
        let elem = ring.add(&a[i], &b[i]);
        c[i] = elem;
    }
    c
}

pub fn poly_sub(ring: &PrimeRing, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    let mut c = vec![ring.zero(); a.len()];
    for i in 0..a.len() {
        let elem = ring.sub(&a[i], &b[i]);
        c[i] = elem;
    }
    c
}

/// Generate the public key `b` for the BGG+ RLWE attribute encoding
/// where `b` is a matrix of ring elements of size `(ell + 1) x m`
/// Using the notation of [DDP+17] each row is the vector b_i
pub fn pub_key_gen(params: &Parameters) -> Vec<Vec<Vec<u64>>> {
    let mut rng = thread_rng();
    let ring = &params.ring;
    let mut b = vec![vec![vec![]; params.m]; params.ell + 1];
    for i in 0..(params.ell + 1) {
        for j in 0..params.m {
            let b_ij = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
            b[i][j] = b_ij;
        }
    }
    b
}

/// Encode the attribute vector `x` into the ciphertext `ct_cin`
///
/// # Arguments
///
/// * `bgg_rlwe`: BGG+ RLWE system
/// * `x`: attribute vector
///
/// # Returns
///
/// * `ct_inner`: ciphertext inner component
/// * `ct_cin`: ciphertext cin component
pub fn encode_attribute_vector(
    bgg_rlwe: &BggRlwe,
    x: &Vec<u64>,
) -> (Vec<Vec<Vec<u64>>>, Vec<Vec<Vec<u64>>>) {
    assert!(x.len() == bgg_rlwe.params.ell + 1);

    let mut rng = thread_rng();
    let ring = &bgg_rlwe.params.ring;
    let s = ring.sample_uniform_vec(ring.ring_size(), &mut rng);
    let mut ct_cin = bgg_rlwe.public_key.clone();
    let mut ct_inner = bgg_rlwe.public_key.clone();

    let mut err_a: Vec<Vec<u64>> = vec![vec![]; bgg_rlwe.params.m];
    let gaussian: NoiseDistribution = Gaussian(3.19).into();

    for i in 0..bgg_rlwe.params.m {
        let mut err = vec![ring.zero(); ring.ring_size()];
        ring.sample_into::<i64>(&mut err, gaussian, rng.clone());
        err_a[i] = err;
    }

    for i in 0..bgg_rlwe.params.ell + 1 {
        // for each i
        // compute the ciphertext error part `err_cin_i`
        let mut err_cin_i = vec![vec![ring.zero(); ring.ring_size()]; bgg_rlwe.params.m];
        for si in 0..bgg_rlwe.params.m {
            for sj in 0..bgg_rlwe.params.m {
                let random_bit = if rng.gen_bool(0.5) { 1 } else { -1 };
                if random_bit == 1 {
                    err_cin_i[si] = poly_add(ring, &err_cin_i[si], &err_a[sj]);
                }
                if random_bit == -1 {
                    err_cin_i[si] = poly_sub(ring, &err_cin_i[si], &err_a[sj]);
                }
            }
        }

        let mut ct_inner_i = bgg_rlwe.public_key[i].clone();

        for j in 0..bgg_rlwe.params.m {
            if x[i] == 1 {
                ct_inner_i[j] = poly_add(ring, &ct_inner_i[j], &bgg_rlwe.params.g[j]);
            }
        }

        let mut ct_cin_i = ct_inner_i.clone();

        for j in 0..bgg_rlwe.params.m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let c = ring.take_poly(&mut scratch);
            ring.poly_mul(c, &ct_inner_i[j], &s, scratch.reborrow());
            ct_cin_i[j] = poly_add(ring, &c.to_vec(), &err_cin_i[j]);
        }

        ct_inner[i] = ct_inner_i;
        ct_cin[i] = ct_cin_i;
    }
    (ct_inner, ct_cin)
}

pub fn m_eval_add(params: &Parameters, bu: &Vec<Vec<u64>>, bv: &Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    let mut bw = vec![vec![]; params.m];
    let ring = &params.ring;
    for i in 0..params.m {
        bw[i] = poly_add(&ring, &bu[i], &bv[i]);
    }
    bw
}

pub fn m_eval_add_x(
    params: &Parameters,
    ct_inner_1: &Vec<Vec<u64>>,
    ct_inner_2: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    let mut ct_inner_1_plus_ct_inner_2 = vec![vec![]; params.m];
    let ring = &params.ring;
    for i in 0..params.m {
        ct_inner_1_plus_ct_inner_2[i] = poly_add(ring, &ct_inner_1[i], &ct_inner_2[i]);
    }
    ct_inner_1_plus_ct_inner_2
}

pub fn m_eval_mul(params: &Parameters, bu: &Vec<Vec<u64>>, bv: &Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    let ring = &params.ring;
    let mut bw = vec![vec![ring.zero(); ring.ring_size()]; params.m];
    // Compute minus_bu by multiplying each coefficient by -1
    let mut minus_bu = vec![vec![ring.zero(); ring.ring_size()]; params.m];
    for i in 0..params.m {
        for j in 0..ring.ring_size() {
            // To get -1 * coefficient in the ring, we subtract the coefficient from 0
            minus_bu[i][j] = ring.sub(&ring.zero(), &bu[i][j]);
        }
    }
    let tau = bit_decompose(params, &minus_bu);

    // Compute bw = bv * TAU
    for i in 0..params.m {
        // For each polynomial in bv
        for h in 0..params.m {
            // For each row of TAU
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let product = ring.take_poly(&mut scratch);

            // Multiply bv[h] by tau[h][i]
            ring.poly_mul(product, &bv[h], &tau[h][i], scratch.reborrow());

            // Add to the result
            bw[i] = poly_add(ring, &bw[i], &product.to_vec());
        }
    }
    bw
}

pub fn m_eval_mul_x(
    params: &Parameters,
    ct_inner_1: &Vec<Vec<u64>>,
    ct_inner_2: &Vec<Vec<u64>>,
    pub_key_1: &Vec<Vec<u64>>,
    x_2: u64,
) -> Vec<Vec<u64>> {
    let mut ct_inner_1_times_ct_inner_2 = vec![vec![]; params.m];
    let ring = &params.ring;

    // First compute x_2 * ct_inner_1
    for i in 0..params.m {
        let mut scratch = ring.allocate_scratch(1, 2, 0);
        let mut scratch = scratch.borrow_mut();
        let ct_inner_1_times_x_2_i = ring.take_poly(&mut scratch);

        // Create polynomial for x_2
        let mut x = vec![ring.zero(); ring.ring_size()];
        x[0] = ring.elem_from(x_2);

        ring.poly_mul(
            ct_inner_1_times_x_2_i,
            &x,
            &ct_inner_1[i],
            scratch.reborrow(),
        );
        ct_inner_1_times_ct_inner_2[i] = ct_inner_1_times_x_2_i.to_vec();
    }

    // Compute τ which is the bit decomposition of pubkey_1
    let mut minus_pub_key_1 = vec![vec![ring.zero(); ring.ring_size()]; params.m];
    for i in 0..params.m {
        for j in 0..ring.ring_size() {
            // To get -1 * coefficient in the ring, we subtract the coefficient from 0
            minus_pub_key_1[i][j] = ring.sub(&ring.zero(), &pub_key_1[i][j]);
        }
    }
    let tau = bit_decompose(params, &minus_pub_key_1);

    // Add τᵀC₂ to the result
    // Note: When accessing tau[i][h] instead of tau[h][i], we're effectively using τᵀ
    for i in 0..params.m {
        for h in 0..params.m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let product = ring.take_poly(&mut scratch);

            // Multiply C₂[h] by τᵀ[i][h] (which is tau[h][i])
            ring.poly_mul(product, &ct_inner_2[h], &tau[h][i], scratch.reborrow());

            // Add to the result
            ct_inner_1_times_ct_inner_2[i] =
                poly_add(ring, &ct_inner_1_times_ct_inner_2[i], &product.to_vec());
        }
    }

    ct_inner_1_times_ct_inner_2
}

pub fn bit_decompose(params: &Parameters, bu: &Vec<Vec<u64>>) -> Vec<Vec<Vec<u64>>> {
    let ring = &params.ring;
    let ring_size = ring.ring_size();
    // Create a matrix of dimension m × m, where each element is a binary polynomial
    let mut tau = vec![vec![vec![ring.zero(); ring_size]; params.m]; params.m];

    // For each row h in the output matrix
    for h in 0..params.m {
        // For each column i in the output matrix
        for i in 0..params.m {
            // For each coefficient j in the polynomial
            for j in 0..ring_size {
                // Get the h-th bit of the j-th coefficient of the i-th polynomial
                let coeff = bu[i][j];
                // Check if the h-th bit is set
                let bit = (coeff >> h) & 1;
                tau[h][i][j] = bit;
            }
        }
    }
    tau
}

#[cfg(test)]
mod tests {
    use phantom_zone_math::prelude::ElemFrom;

    use super::*;
    use crate::BggRlwe;

    #[test]
    fn test_matrix_encoding_homomorphism_add_gate() {
        let bgg_rlwe = BggRlwe::new(12, 51, 4);
        let mut rng = thread_rng();
        let mut x = (0..bgg_rlwe.params.ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let (ct_inner, _) = encode_attribute_vector(&bgg_rlwe, &x);

        // Perform plus gate of b[1] and b[2]
        let h_1_plus_2 = m_eval_add(
            &bgg_rlwe.params,
            &bgg_rlwe.public_key[1],
            &bgg_rlwe.public_key[2],
        );

        // Perform plus gate of ct_inner[1] and ct_inner[2]
        let h_1_plus_2_x = m_eval_add_x(&bgg_rlwe.params, &ct_inner[1], &ct_inner[2]);

        // Verify homomorphism of plus gate
        let lhs = h_1_plus_2_x;

        let ring = &bgg_rlwe.params.ring;
        let mut rhs = h_1_plus_2.clone();
        // Create fx such that the constant term (first coefficient) is to x[1] + x[2]
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] + x[2]);

        for i in 0..bgg_rlwe.params.m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &bgg_rlwe.params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..bgg_rlwe.params.m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }

    #[test]
    fn test_matrix_encoding_homomorphism_mul_gate() {
        let bgg_rlwe = BggRlwe::new(12, 51, 4);
        let mut rng = thread_rng();
        let mut x = (0..bgg_rlwe.params.ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let (ct_inner, _) = encode_attribute_vector(&bgg_rlwe, &x);

        // Perform multiplication gate of b[1] and b[2]
        let h_1_times_2 = m_eval_mul(
            &bgg_rlwe.params,
            &bgg_rlwe.public_key[1],
            &bgg_rlwe.public_key[2],
        );

        // Perform multiplication gate of ct_inner[1] and ct_inner[2]
        let h_1_times_2_x = m_eval_mul_x(
            &bgg_rlwe.params,
            &ct_inner[1],
            &ct_inner[2],
            &bgg_rlwe.public_key[1],
            x[2], // Pass x[2] as the scalar multiplier
        );

        // Verify homomorphism of multiplication gate
        let lhs = h_1_times_2_x;

        let ring = &bgg_rlwe.params.ring;
        let mut rhs = h_1_times_2.clone();
        // Create fx such that the constant term (first coefficient) is x[1] * x[2]
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] * x[2]);

        for i in 0..bgg_rlwe.params.m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &bgg_rlwe.params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }

        for i in 0..bgg_rlwe.params.m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }

    #[test]
    fn test_bit_decompose() {
        let bgg_rlwe = BggRlwe::new(12, 51, 4);
        let bu = bgg_rlwe.public_key[1].clone();
        let ring = &bgg_rlwe.params.ring;
        let tau = bit_decompose(&bgg_rlwe.params, &bu);

        // Reconstruct the original input by multiplying tau with G
        let mut reconstructed = vec![vec![ring.zero(); ring.ring_size()]; bgg_rlwe.params.m];

        // For each column i of the output
        for i in 0..bgg_rlwe.params.m {
            // For each row h of tau
            for h in 0..bgg_rlwe.params.m {
                // Multiply tau[h][i] by g[h] and add to the result
                let mut scratch = ring.allocate_scratch(1, 2, 0);
                let mut scratch = scratch.borrow_mut();
                let product = ring.take_poly(&mut scratch);
                ring.poly_mul(
                    product,
                    &tau[h][i],
                    &bgg_rlwe.params.g[h],
                    scratch.reborrow(),
                );
                reconstructed[i] = poly_add(ring, &reconstructed[i], &product.to_vec());
            }
        }

        for i in 0..bgg_rlwe.params.m {
            assert_eq!(bu[i], reconstructed[i]);
        }
    }
}
