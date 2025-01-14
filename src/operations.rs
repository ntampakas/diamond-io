use phantom_zone_crypto::util::distribution::NoiseDistribution;
use phantom_zone_math::{
    prelude::{Gaussian, ModulusOps, Sampler},
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
        let mut err_cin_i: Vec<Vec<u64>> =
            vec![vec![ring.zero(); ring.ring_size()]; bgg_rlwe.params.m];
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
    let mut ct_inner_1_plus_2 = vec![vec![]; params.m];
    let ring = &params.ring;
    for i in 0..params.m {
        ct_inner_1_plus_2[i] = poly_add(ring, &ct_inner_1[i], &ct_inner_2[i]);
    }
    ct_inner_1_plus_2
}

#[cfg(test)]
mod tests {
    use phantom_zone_math::prelude::ElemFrom;

    use super::*;
    use crate::BggRlwe;

    #[test]
    fn test_matrix_encoding_homomorphism() {
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
}
