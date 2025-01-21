use phantom_zone_math::{prelude::ModulusOps, ring::RingOps};

use crate::{
    operations::{bit_decompose, gen_identity_matrix_to_scalar, poly_add},
    parameters::Parameters,
    utils::empty_vector_ring,
};

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

pub fn m_eval_mul(
    params: &Parameters,
    b_left: &Vec<Vec<u64>>,
    b_right: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    let ring = params.ring();
    let m = *params.m();
    let mut out = empty_vector_ring(ring, m);

    // Compute minus_b_left by multiplying each coefficient by -1
    let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            // To get -1 * coefficient in the ring, we subtract the coefficient from 0
            minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
        }
    }

    let tau = bit_decompose(params, &minus_b_left);

    // Compute out = b_right * TAU
    for i in 0..m {
        for h in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let product = ring.take_poly(&mut scratch);

            // Multiply b_right[h] by tau[h][i]
            ring.poly_mul(product, &b_right[h], &tau[h][i], scratch.reborrow());

            out[i] = poly_add(ring, &out[i], &product.to_vec());
        }
    }
    out
}

pub fn m_eval_add_x(
    params: &Parameters,
    b_left: &Vec<Vec<u64>>,
    x_left: &u64,
    b_right: &Vec<Vec<u64>>,
    x_right: &u64,
) -> (Vec<Vec<Vec<u64>>>, Vec<Vec<Vec<u64>>>) {
    let ring = params.ring();
    let m = *params.m();
    let out_left = gen_identity_matrix_to_scalar(ring, m, 1);
    let out_right = gen_identity_matrix_to_scalar(ring, m, 1);

    (out_left, out_right)
}

pub fn m_eval_mul_x(
    params: &Parameters,
    b_left: &Vec<Vec<u64>>,
    x_left: &u64,
    b_right: &Vec<Vec<u64>>,
    x_right: &u64,
) -> (Vec<Vec<Vec<u64>>>, Vec<Vec<Vec<u64>>>) {
    let ring = params.ring();
    let m = *params.m();

    // First matrix: Identity matrix scaled by x_right
    let out_left = gen_identity_matrix_to_scalar(ring, m, *x_right);

    // Second matrix: Tau(-b_left)
    // First compute -b_left
    let mut minus_b_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            minus_b_left[i][j] = ring.sub(&ring.zero(), &b_left[i][j]);
        }
    }

    // Compute tau of -b_left
    let out_right = bit_decompose(params, &minus_b_left);

    (out_left, out_right)
}

#[cfg(test)]
mod tests {
    use crate::{
        ciphertext::Ciphertext,
        eval::{m_eval_add, m_eval_add_x, m_eval_mul, m_eval_mul_x},
        operations::{mat_vert_concat, poly_add, vec_mat_mul},
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
        let ct_inner = ciphertext.inner();

        // Perform add gate of b[1] and b[2]
        let b_1_plus_2 = m_eval_add(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

        // Perform add gate of ct_inner[1] and ct_inner[2]
        let (scalar_left, scalar_right) = m_eval_add_x(
            pub_key.params(),
            &pub_key.b()[1],
            &x[1],
            &pub_key.b()[2],
            &x[2],
        );

        // Define h_1_plus_2_x as the vertical concatenation of scalar_left and scalar_right
        let h_1_plus_2_x = mat_vert_concat(ring, &scalar_left, &scalar_right);

        // Verify homomorphism of add gate such that (ct_inner[1] | ct_inner[2]) * h_1_plus_2_x = b_1_plus_2 + (x1+x2)G
        let concat_vec = [ct_inner[1].clone(), ct_inner[2].clone()].concat();
        let lhs = vec_mat_mul(ring, &concat_vec, &h_1_plus_2_x);

        let mut rhs = b_1_plus_2;
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
        let ct_inner = ciphertext.inner();

        // Perform mul gate of b[1] and b[2]
        let b_1_times_2 = m_eval_mul(pub_key.params(), &pub_key.b()[1], &pub_key.b()[2]);

        // Perform mul gate of ct_inner[1] and ct_inner[2]
        let (scalar_left, scalar_right) = m_eval_mul_x(
            pub_key.params(),
            &pub_key.b()[1],
            &x[1],
            &pub_key.b()[2],
            &x[2],
        );

        // Define h_1_times_2_x as the vertical concatenation of scalar_left and scalar_right
        let h_1_times_2_x = mat_vert_concat(ring, &scalar_left, &scalar_right);

        // Verify homomorphism of mul gate such that (ct_inner[1] | ct_inner[2]) * h_1_times_2_x = b_1_times_2 + (x1*x2)G
        let concat_vec = [ct_inner[1].clone(), ct_inner[2].clone()].concat();
        let lhs = vec_mat_mul(ring, &concat_vec, &h_1_times_2_x);

        let mut rhs = b_1_times_2.clone();
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
