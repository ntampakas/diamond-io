mod ciphertext;
mod operations;
mod parameters;
mod pub_key;
mod utils;

pub use ciphertext::Ciphertext;
pub use operations::{bit_decompose, poly_add, poly_sub};
pub use parameters::Parameters;
pub use pub_key::PublicKey;
pub use utils::{print_matrix_ring, print_vector, print_vector_ring};

#[cfg(test)]
mod tests {
    use crate::{operations::vec_mat_mul, poly_add, Ciphertext, Parameters, PublicKey};
    use phantom_zone_math::{
        prelude::{ElemFrom, ModulusOps},
        ring::RingOps,
    };
    use rand::{thread_rng, Rng};

    #[test]
    fn test_matrix_encoding_homomorphism_add_gate() {
        let params = Parameters::new(12, 51, 4);
        let mut pub_key = PublicKey::new(&params);
        let mut rng = thread_rng();
        let ring = &params.ring;
        let m = params.m;
        let mut x = (0..params.ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let ciphertext = Ciphertext::new(&pub_key, &params, &x);
        let ct_inner = ciphertext.inner();

        // Perform add gate of b[1] and b[2]
        pub_key.add_gate(1, 2);
        let b_1_plus_2 = pub_key.b()[params.ell + 1].clone();

        // Perform add gate of ct_inner[1] and ct_inner[2]
        let h_1_plus_2_x = ciphertext.add_gate(1, 2);

        // Verify homomorphism of add gate such that (ct_inner[1] | ct_inner[2]) * h_1_plus_2_x = b_1_plus_2 + (x1+x2)G
        let concat_vec = [ct_inner[1].clone(), ct_inner[2].clone()].concat();
        let lhs = vec_mat_mul(ring, concat_vec, h_1_plus_2_x);

        let mut rhs = b_1_plus_2;
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] + x[2]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }

    #[test]
    fn test_matrix_encoding_homomorphism_mul_gate() {
        let params = Parameters::new(12, 51, 4);
        let mut pub_key = PublicKey::new(&params);
        let mut rng = thread_rng();
        let ring = &params.ring;
        let m = params.m;
        let mut x = (0..params.ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let ciphertext = Ciphertext::new(&pub_key, &params, &x);
        let ct_inner = ciphertext.inner();

        // Perform mul gate of b[1] and b[2]
        pub_key.mul_gate(1, 2);
        let b_1_times_2 = pub_key.b()[params.ell + 1].clone();

        // Perform add gate of ct_inner[1] and ct_inner[2]
        let h_1_times_2_x = ciphertext.mul_gate(1, 2);

        // Verify homomorphism of mul gate such that (ct_inner[1] | ct_inner[2]) * h_1_times_2_x = b_1_times_2 + (x1*x2)G
        let concat_vec = [ct_inner[1].clone(), ct_inner[2].clone()].concat();
        let lhs = vec_mat_mul(ring, concat_vec, h_1_times_2_x);

        let mut rhs = b_1_times_2.clone();
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] * x[2]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }

    #[test]
    fn test_matrix_encoding_homomorphism_circuit() {
        let params = Parameters::new(12, 51, 7);
        let mut pub_key = PublicKey::new(&params);
        let mut rng = thread_rng();
        let ring = &params.ring;
        let m = params.m;
        let mut x = (0..params.ell + 1)
            .map(|_| rng.gen_range(0..2))
            .collect::<Vec<_>>();
        x[0] = 1; // The actual attribute vector is x[1..], the value set to the index 0 is just for easier arithmetic during encoding

        let ciphertext = Ciphertext::new(&pub_key, &params, &x);
        let ct_inner = ciphertext.inner();

        // Circuit: (x1 + x2) * (x3)
        // Perform add gate of b[1] and b[2]
        pub_key.add_gate(1, 2);
        let b_gate_first = pub_key.b()[params.ell + 1].clone();

        // Perform add gate of ct_inner[1] and ct_inner[2]
        let h_gate_first_x = ciphertext.add_gate(1, 2);

        // Verify homomorphism of add gate such that (ct_inner[1] | ct_inner[2]) * h_gate_first_x = b_gate_first + (x1+x2)G
        let concat_vec = [ct_inner[1].clone(), ct_inner[2].clone()].concat();
        let lhs = vec_mat_mul(ring, concat_vec, h_gate_first_x);

        let mut rhs = b_gate_first;
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] + x[2]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }

        // Perform add gate of b_gate_first and b[3]
        pub_key.add_gate(params.ell + 1, 3);
        let b_gate_second = pub_key.b()[params.ell + 1].clone();

        // Perform add gate of ct_inner[params.ell + 1] and ct_inner[3]
        let h_gate_second_x = ciphertext.add_gate(params.ell + 1, 3);

        // Verify homomorphism of mul gate such that (b_gate_first + (x1+x2)G | ct_inner[3]) * h_gate_second_x = b_gate_second + ((x1+x2+x3)G
        let concat_vec = [rhs.clone(), ct_inner[3].clone()].concat();
        let lhs = vec_mat_mul(ring, concat_vec, h_gate_second_x);

        let mut rhs = b_gate_second;
        let mut fx = vec![ring.zero(); ring.ring_size()];
        fx[0] = ring.elem_from(x[1] + x[2] + x[3]);

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let gi_times_fx = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx, &params.g[i], &fx, scratch.reborrow());
            let gi_times_fx_vec = gi_times_fx.to_vec();
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx_vec);
        }
        for i in 0..m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }
}
