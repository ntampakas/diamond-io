use phantom_zone_math::{prelude::ModulusOps, ring::RingOps};

use crate::{
    operations::{
        bit_decompose, gen_identity_matrix_to_scalar, poly_add, vec_mat_mul, vec_vec_add,
    },
    parameters::Parameters,
    utils::empty_vector_ring,
};

/// `out = b_left + b_right`
pub fn m_eval_add(
    params: &Parameters,
    b_left: &Vec<Vec<u64>>,
    b_right: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    let ring = params.ring();
    let m = *params.m();
    let mut out = empty_vector_ring(ring, m);
    for i in 0..m {
        out[i] = poly_add(ring, &b_left[i], &b_right[i]);
    }

    out
}

/// `out = b_right * tau(-b_left)`
pub fn m_eval_mul(
    params: &Parameters,
    b_left: &Vec<Vec<u64>>,
    b_right: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    let ring = params.ring();
    let m = *params.m();
    let mut out = empty_vector_ring(ring, m);

    let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            // To get -1 * coefficient in the ring, we subtract the coefficient from 0
            minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
        }
    }

    let tau = bit_decompose(params, &minus_b_left);

    for i in 0..m {
        for h in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let product = ring.take_poly(&mut scratch);

            ring.poly_mul(product, &b_right[h], &tau[h][i], scratch.reborrow());

            out[i] = poly_add(ring, &out[i], &product.to_vec());
        }
    }

    out
}

/// `out = ct_left + ct_right`
pub fn m_eval_add_x(
    params: &Parameters,
    ct_left: &Vec<Vec<u64>>,
    ct_right: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    let ring = params.ring();
    vec_vec_add(ring, ct_left, ct_right)
}

/// `out = ct_left * identity_matrix_scaled_by_x_right + ct_right * tau(-b_left)`
pub fn m_eval_mul_x(
    params: &Parameters,
    ct_left: &Vec<Vec<u64>>,
    b_left: &Vec<Vec<u64>>,
    ct_right: &Vec<Vec<u64>>,
    x_right: &u64,
) -> Vec<Vec<u64>> {
    let ring = params.ring();
    let m = *params.m();

    let part_left = gen_identity_matrix_to_scalar(ring, m, *x_right);

    let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
        }
    }

    let part_right = bit_decompose(params, &minus_b_left);

    let out_left = vec_mat_mul(ring, ct_left, &part_left); // TODO: make it more efficient

    let out_right = vec_mat_mul(ring, ct_right, &part_right);
    let out = vec_vec_add(ring, &out_left, &out_right);

    out
}

#[cfg(test)]
mod tests {
    use crate::{
        ciphertext::Ciphertext,
        eval::{m_eval_add, m_eval_add_x, m_eval_mul, m_eval_mul_x},
        operations::poly_add,
        parameters::Parameters,
        pub_key::PublicKey,
    };
    use phantom_zone_math::{
        prelude::{ElemFrom, ModulusOps},
        ring::RingOps,
    };
    use rand::{thread_rng, Rng};

    #[test]
    fn test_matrix_encoding_homomorphism_add_gate() {
        let params = Parameters::new(12, 51, 7);
        let pub_key = PublicKey::new(params);
        let mut rng = thread_rng();
        let ring = pub_key.params().ring();
        let m = *pub_key.params().m();
        let ell = *pub_key.params().ell();
        let g = pub_key.params().g();
        let mut x = (0..ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let ciphertext = Ciphertext::new(&pub_key, &x);

        // Perform add gate of b[1] and b[2]
        let b_out = m_eval_add(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

        // Perform add gate of ct_inner[1] and ct_inner[2]
        let ct_out = m_eval_add_x(
            pub_key.params(),
            &ciphertext.inner()[1],
            &ciphertext.inner()[2],
        );

        let lhs = ct_out.clone();

        // rhs = b_out + g * (x1 + x2)
        let mut rhs = b_out.clone();
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] + x[2]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }

    #[test]
    fn test_matrix_encoding_homomorphism_mul_gate() {
        let params = Parameters::new(12, 51, 7);
        let pub_key = PublicKey::new(params);
        let mut rng = thread_rng();
        let ring = pub_key.params().ring();
        let m = *pub_key.params().m();
        let ell = *pub_key.params().ell();
        let g = pub_key.params().g();
        let mut x = (0..ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let ciphertext = Ciphertext::new(&pub_key, &x);

        // Perform mul gate of b[1] and b[2]
        let b_out = m_eval_mul(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

        // Perform mul gate of ct_inner[1] and ct_inner[2]
        let ct_out = m_eval_mul_x(
            pub_key.params(),
            &ciphertext.inner()[1],
            &pub_key.b()[1],
            &ciphertext.inner()[2],
            &x[2],
        );

        let lhs = ct_out.clone();

        // rhs = b_out + g * (x1 * x2)
        let mut rhs = b_out.clone();
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] * x[2]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }
}
