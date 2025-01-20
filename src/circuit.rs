use phantom_zone_math::{
    prelude::{ElemFrom, ModulusOps},
    ring::RingOps,
};

use crate::{
    eval::{m_eval_add, m_eval_add_x, m_eval_mul, m_eval_mul_x},
    operations::{mat_mat_add, mat_mat_mul},
    parameters::Parameters,
    pub_key::PublicKey,
};

#[derive(Clone, Debug)]
pub struct CircuitX {
    b_gates: Vec<Vec<Vec<u64>>>,
    x_gates: Vec<u64>,
    h_gates: Vec<Vec<Vec<Vec<u64>>>>,
}

impl CircuitX {
    pub fn new(pub_key: &PublicKey, x: &Vec<u64>) -> Self {
        let params = pub_key.params();
        let ring = params.ring();
        let m = *params.m();
        let ell = *params.ell();

        let b_gates = pub_key.b().clone();
        let x_gates = x.clone();
        let mut h_gates = Vec::with_capacity(ell + 1);

        for i in 0..ell + 1 {
            let mut h_gate = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; (ell + 1) * m];
            for j in 0..m {
                h_gate[i * m + j][j][0] = ring.elem_from(1u64);
            }

            h_gates.push(h_gate);
        }

        Self {
            b_gates,
            x_gates,
            h_gates,
        }
    }

    pub fn b_gates(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.b_gates
    }

    pub fn x_gates(&self) -> &Vec<u64> {
        &self.x_gates
    }

    pub fn h_gates(&self) -> &Vec<Vec<Vec<Vec<u64>>>> {
        &self.h_gates
    }

    pub fn add_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
        let ring = params.ring();

        // Calculate b for the new gate
        let b_gate = m_eval_add(
            params,
            &self.b_gates()[idx_left],
            &self.b_gates()[idx_right],
        );
        self.b_gates.push(b_gate);

        // Calculate x for the new gate
        let x_gate = self.x_gates()[idx_left] + self.x_gates()[idx_right];
        self.x_gates.push(x_gate);

        // Calculate h for the new gate
        let (scalar_left, scalar_right) = m_eval_add_x(
            params,
            &self.b_gates()[idx_left],
            &self.x_gates()[idx_left],
            &self.b_gates()[idx_right],
            &self.x_gates()[idx_right],
        );

        let h_input_left = &self.h_gates()[idx_left];
        let h_input_right = &self.h_gates()[idx_right];

        let mat_a = mat_mat_mul(ring, h_input_left, &scalar_left);
        let mat_b = mat_mat_mul(ring, h_input_right, &scalar_right);

        let h_gate = mat_mat_add(ring, &mat_a, &mat_b);

        self.h_gates.push(h_gate);
        self.h_gates.len() - 1
    }

    pub fn mul_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
        let ring = params.ring();

        // Calculate b for the new gate
        let b_gate = m_eval_mul(
            params,
            &self.b_gates()[idx_left],
            &self.b_gates()[idx_right],
        );
        self.b_gates.push(b_gate);

        // Calculate x for the new gate
        let x_gate = self.x_gates()[idx_left] * self.x_gates()[idx_right];
        self.x_gates.push(x_gate);

        // Calculate h for the new gate
        let (scalar_left, scalar_right) = m_eval_mul_x(
            params,
            &self.b_gates()[idx_left],
            &self.x_gates()[idx_left],
            &self.b_gates()[idx_right],
            &self.x_gates()[idx_right],
        );

        let h_input_left = &self.h_gates()[idx_left];
        let h_input_right = &self.h_gates()[idx_right];

        let mat_a = mat_mat_mul(ring, h_input_left, &scalar_left);
        let mat_b = mat_mat_mul(ring, h_input_right, &scalar_right);

        let h_gate = mat_mat_add(ring, &mat_a, &mat_b);

        self.h_gates.push(h_gate);
        self.h_gates.len() - 1
    }
}

#[derive(Clone, Debug)]
pub struct Circuit {
    b_gates: Vec<Vec<Vec<u64>>>,
}

impl Circuit {
    pub fn new(pub_key: &PublicKey) -> Self {
        let b_gates = pub_key.b().clone();
        Self { b_gates }
    }

    pub fn b_gates(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.b_gates
    }

    pub fn add_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
        // Calculate b for the new gate
        let b_gate = m_eval_add(
            params,
            &self.b_gates()[idx_left],
            &self.b_gates()[idx_right],
        );
        self.b_gates.push(b_gate);
        self.b_gates.len() - 1
    }

    pub fn mul_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
        let b_gate = m_eval_mul(
            params,
            &self.b_gates()[idx_left],
            &self.b_gates()[idx_right],
        );
        self.b_gates.push(b_gate);
        self.b_gates.len() - 1
    }
}

#[cfg(test)]
mod tests {
    use crate::ciphertext::Ciphertext;
    use crate::circuit::{Circuit, CircuitX};
    use crate::operations::{mat_horiz_concat, poly_add, vec_horiz_concat, vec_mat_mul};
    use crate::parameters::Parameters;
    use crate::pub_key::PublicKey;
    use phantom_zone_math::prelude::{ElemFrom, ModulusOps};
    use phantom_zone_math::ring::RingOps;
    use rand::{thread_rng, Rng};

    #[test]
    fn test_circuit_single_output() {
        let params = Parameters::new(3, 12, 7);
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
        let ct_inner_concat = ciphertext.inner_concat();

        let mut circuit_x = CircuitX::new(&pub_key, &x);
        let mut circuit = Circuit::new(&pub_key);

        let random_bit_gate_4 = rng.gen_range(0..2);
        let random_bit_gate_5 = rng.gen_range(0..2);
        let gate_idx_4: usize;
        let gate_idx_5: usize;

        if random_bit_gate_4 == 1 {
            // x4 = x1 * x2
            gate_idx_4 = circuit.mul_gate(pub_key.params(), 1, 2);
            circuit_x.mul_gate(pub_key.params(), 1, 2);
        } else {
            // x4 = x1 + x2
            gate_idx_4 = circuit.add_gate(pub_key.params(), 1, 2);
            circuit_x.add_gate(pub_key.params(), 1, 2);
        }

        if random_bit_gate_5 == 1 {
            // x5 = x4 * x3
            gate_idx_5 = circuit.mul_gate(pub_key.params(), gate_idx_4, 3);
            circuit_x.mul_gate(pub_key.params(), gate_idx_4, 3);
        } else {
            // x5 = x4 + x3
            gate_idx_5 = circuit.add_gate(pub_key.params(), gate_idx_4, 3);
            circuit_x.add_gate(pub_key.params(), gate_idx_4, 3);
        }

        // verify homomorphism such that (ct_inner[0] | ct_inner[1] | ct_inner[2] | ... | ct_inner[ell]) * h[gate_idx_5] = b[gate_idx_5] + f(x)G
        let lhs = vec_mat_mul(ring, &ct_inner_concat, &circuit_x.h_gates()[gate_idx_5]);

        let mut rhs = circuit.b_gates()[gate_idx_5].clone();
        let mut fx = vec![ring.zero(); ring.ring_size()];

        if random_bit_gate_4 == 1 && random_bit_gate_5 == 1 {
            fx[0] = ring.elem_from((x[1] * x[2]) * x[3]);
        } else if random_bit_gate_4 == 1 && random_bit_gate_5 == 0 {
            fx[0] = ring.elem_from((x[1] * x[2]) + x[3]);
        } else if random_bit_gate_4 == 0 && random_bit_gate_5 == 1 {
            fx[0] = ring.elem_from((x[1] + x[2]) * x[3]);
        } else {
            fx[0] = ring.elem_from((x[1] + x[2]) + x[3]);
        }

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
    fn test_circuit_add_multi_output() {
        let params = Parameters::new(3, 12, 7);
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
        let ct_inner_concat = ciphertext.inner_concat();

        let mut circuit_x = CircuitX::new(&pub_key, &x);
        let mut circuit = Circuit::new(&pub_key);

        let gate_idx_4 = circuit.add_gate(pub_key.params(), 1, 2); // x1 + x2
        let _ = circuit.mul_gate(pub_key.params(), gate_idx_4, 3); // (x1 + x2) * x3
        let _ = circuit.mul_gate(pub_key.params(), gate_idx_4, 4); // (x1 + x2) * x4

        let gate_idx_4 = circuit_x.add_gate(pub_key.params(), 1, 2); // x1 + x2
        let gate_idx_5 = circuit_x.mul_gate(pub_key.params(), gate_idx_4, 3); // (x1 + x2) * x3
        let gate_idx_6 = circuit_x.mul_gate(pub_key.params(), gate_idx_4, 4); // (x1 + x2) * x4

        // verify homomorphism such that (ct_inner[0] | ct_inner[1] | ct_inner[2] | ... | ct_inner[ell]) * (h[gate_idx_5] | h[gate_idx_6]) = [b[gate_idx_5] | b[gate_idx_6]] + [((x1 + x2) * x3)G | ((x1 + x2) * x4)G]

        // horizontally concatenate the two matrices h[gate_idx_5] and h[gate_idx_6]
        let concatenated_h_gates = mat_horiz_concat(
            ring,
            &circuit_x.h_gates()[gate_idx_5],
            &circuit_x.h_gates()[gate_idx_6],
        );

        let lhs = vec_mat_mul(ring, &ct_inner_concat, &concatenated_h_gates);

        // define rhs as the concatenation of b[gate_idx_5] and b[gate_idx_6]
        let mut rhs = vec_horiz_concat(
            &circuit.b_gates()[gate_idx_5],
            &circuit.b_gates()[gate_idx_6],
        );

        // Add (x1 + x2 + x3)G to the left half and (x1 + x2 + x4)G to the right half
        let mut fx1 = vec![ring.zero(); ring.ring_size()];
        let mut fx2 = vec![ring.zero(); ring.ring_size()];
        fx1[0] = ring.elem_from((x[1] + x[2]) * x[3]); // (x1 + x2) * x3
        fx2[0] = ring.elem_from((x[1] + x[2]) * x[4]); // (x1 + x2) * x4

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();

            // For the left half (gate_idx_5)
            let gi_times_fx1 = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx1, &g[i], &fx1, scratch.reborrow());
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx1.to_vec());

            // For the right half (gate_idx_6)
            let gi_times_fx2 = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx2, &g[i], &fx2, scratch.reborrow());
            rhs[i + m] = poly_add(ring, &rhs[i + m], &gi_times_fx2.to_vec());
        }

        for i in 0..2 * m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }
}
