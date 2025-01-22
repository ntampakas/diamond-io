use crate::{
    ciphertext::Ciphertext,
    eval::{m_eval_add, m_eval_add_x, m_eval_mul, m_eval_mul_x},
    parameters::Parameters,
    pub_key::PublicKey,
};

#[derive(Clone, Debug)]
pub struct CircuitX {
    b_gates: Vec<Vec<Vec<u64>>>,
    x_gates: Vec<u64>,
    ct_gates: Vec<Vec<Vec<u64>>>,
}

impl CircuitX {
    pub fn new(pub_key: &PublicKey, x: &Vec<u64>, ciphertext: &Ciphertext) -> Self {
        let params = pub_key.params();
        let ell = *params.ell();

        let b_gates = pub_key.b().clone();
        let x_gates = x.clone();
        let mut ct_gates = Vec::with_capacity(ell + 1);

        for i in 0..ell + 1 {
            let ct_gate = ciphertext.inner()[i].clone(); // TODO: eventually we will need to take the full ciphertext here
            ct_gates.push(ct_gate);
        }

        Self {
            b_gates,
            x_gates,
            ct_gates,
        }
    }

    pub fn b_gates(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.b_gates
    }

    pub fn x_gates(&self) -> &Vec<u64> {
        &self.x_gates
    }

    pub fn ct_gates(&self) -> &Vec<Vec<Vec<u64>>> {
        &self.ct_gates
    }

    pub fn add_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
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

        // Calculate ct for the new gate
        let ct_gate = m_eval_add_x(
            params,
            &self.ct_gates()[idx_left],
            &self.ct_gates()[idx_right],
        );
        self.ct_gates.push(ct_gate);
        self.ct_gates.len() - 1
    }

    pub fn mul_gate(&mut self, params: &Parameters, idx_left: usize, idx_right: usize) -> usize {
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

        // Calculate ct for the new gate
        let ct_gate = m_eval_mul_x(
            params,
            &self.ct_gates()[idx_left],
            &self.b_gates()[idx_left],
            &self.ct_gates()[idx_right],
            &self.x_gates()[idx_right],
        );
        self.ct_gates.push(ct_gate);
        self.ct_gates.len() - 1
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
        // Calculate b for the new gate
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
    use crate::operations::{poly_add, vec_horiz_concat};
    use crate::parameters::Parameters;
    use crate::pub_key::PublicKey;
    use phantom_zone_math::prelude::{ElemFrom, ModulusOps};
    use phantom_zone_math::ring::RingOps;
    use rand::{thread_rng, Rng};

    #[test]
    fn test_circuit_single_output() {
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

        let mut circuit_x = CircuitX::new(&pub_key, &x, &ciphertext);
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

        // verify homomorphism such that ct_out = b_out + f(x)G
        let lhs = circuit_x.ct_gates()[gate_idx_5].clone();

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
    fn test_circuit_multi_output() {
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

        let mut circuit_x = CircuitX::new(&pub_key, &x, &ciphertext);
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
            // x5 = x2 * x3
            gate_idx_5 = circuit.mul_gate(pub_key.params(), 2, 3);
            circuit_x.mul_gate(pub_key.params(), 2, 3);
        } else {
            // x5 = x2 + x3
            gate_idx_5 = circuit.add_gate(pub_key.params(), 2, 3);
            circuit_x.add_gate(pub_key.params(), 2, 3);
        }

        // verify homomorphism such that ct_out = b_out + f(x)G
        let lhs = vec_horiz_concat(
            &circuit_x.ct_gates()[gate_idx_4],
            &circuit_x.ct_gates()[gate_idx_5],
        );

        let mut rhs = vec_horiz_concat(
            &circuit.b_gates()[gate_idx_4],
            &circuit.b_gates()[gate_idx_5],
        );

        let mut fx1 = vec![ring.zero(); ring.ring_size()];
        let mut fx2 = vec![ring.zero(); ring.ring_size()];

        if random_bit_gate_4 == 1 {
            fx1[0] = ring.elem_from(x[1] * x[2]); // x1 * x2
        } else {
            fx1[0] = ring.elem_from(x[1] + x[2]); // x1 + x2
        }

        if random_bit_gate_5 == 1 {
            fx2[0] = ring.elem_from(x[2] * x[3]); // x2 * x3
        } else {
            fx2[0] = ring.elem_from(x[2] + x[3]); // x2 + x3
        }

        for i in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();

            // For the left half g*fx1
            let gi_times_fx1 = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx1, &g[i], &fx1, scratch.reborrow());
            rhs[i] = poly_add(ring, &rhs[i], &gi_times_fx1.to_vec());

            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();

            // For the right half g*fx2
            let gi_times_fx2 = ring.take_poly(&mut scratch);
            ring.poly_mul(gi_times_fx2, &g[i], &fx2, scratch.reborrow());
            rhs[i + m] = poly_add(ring, &rhs[i + m], &gi_times_fx2.to_vec());
        }

        for i in 0..2 * m {
            assert_eq!(lhs[i], rhs[i]);
        }
    }
}
