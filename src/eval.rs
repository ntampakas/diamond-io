use phantom_zone_math::{
    prelude::{ElemFrom, ModulusOps},
    ring::RingOps,
};

use crate::{
    operations::{bit_decompose, poly_add},
    pub_key::PublicKey,
};

pub fn m_eval_add(pub_key: &PublicKey, idx_left: usize, idx_right: usize) -> Vec<Vec<u64>> {
    let ring = pub_key.params().ring();
    let m = *pub_key.params().m();
    let mut out = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        out[i] = poly_add(&ring, &pub_key.b()[idx_left][i], &pub_key.b()[idx_right][i]);
    }
    out
}

pub fn m_eval_mul(pub_key: &PublicKey, idx_left: usize, idx_right: usize) -> Vec<Vec<u64>> {
    let ring = pub_key.params().ring();
    let m = *pub_key.params().m();
    let mut out = vec![vec![ring.zero(); ring.ring_size()]; m];

    // Compute minus_b_idx_left by multiplying each coefficient by -1
    let mut minus_b_idx_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            // To get -1 * coefficient in the ring, we subtract the coefficient from 0
            minus_b_idx_left[i][j] = ring.sub(&ring.zero(), &pub_key.b()[idx_left][i][j]);
        }
    }

    let tau = bit_decompose(pub_key.params(), &minus_b_idx_left);

    // Compute out = b_idx_right * TAU
    for i in 0..m {
        for h in 0..m {
            let mut scratch = ring.allocate_scratch(1, 2, 0);
            let mut scratch = scratch.borrow_mut();
            let product = ring.take_poly(&mut scratch);

            // Multiply b_idx_right[h] by tau[h][i]
            ring.poly_mul(
                product,
                &pub_key.b()[idx_right][h],
                &tau[h][i],
                scratch.reborrow(),
            );

            out[i] = poly_add(ring, &out[i], &product.to_vec());
        }
    }
    out
}

pub fn m_eval_add_x(
    pub_key: &PublicKey,
    x: &Vec<u64>,
    idx_left: usize,
    idx_right: usize,
) -> Vec<Vec<Vec<u64>>> {
    let params = pub_key.params();
    let ring = params.ring();
    let m = *params.m();
    let mut out = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; 2 * m];

    // Fill both identity matrices at once
    for i in 0..m {
        // Set constant polynomial 1 on diagonals of both identity matrices
        out[i][i][0] = ring.elem_from(1u64);
        out[i + m][i][0] = ring.elem_from(1u64);
    }

    out
}

pub fn m_eval_mul_x(
    pub_key: &PublicKey,
    x: &Vec<u64>,
    idx_left: usize,
    idx_right: usize,
) -> Vec<Vec<Vec<u64>>> {
    let params = pub_key.params();
    let ring = params.ring();
    let m = *params.m();

    let mut out = vec![vec![vec![ring.zero(); ring.ring_size()]; m]; 2 * m];

    // First matrix: Identity matrix scaled by x_idx_right
    let mut poly = vec![ring.zero(); ring.ring_size()];
    let x_idx_right = x[idx_right];
    poly[0] = ring.elem_from(x_idx_right);
    for i in 0..m {
        out[i][i] = poly.clone();
    }

    // Second matrix: Tau(-b_idx_left)
    // First compute -b_idx_left
    let b_idx_left = &pub_key.b()[idx_left];
    let mut minus_b_idx_left = vec![vec![ring.zero(); ring.ring_size()]; m];
    for i in 0..m {
        for j in 0..ring.ring_size() {
            minus_b_idx_left[i][j] = ring.sub(&ring.zero(), &b_idx_left[i][j]);
        }
    }

    // Compute tau of -b_idx_left
    let tau = bit_decompose(params, &minus_b_idx_left);

    // Copy tau into the bottom half of out
    for h in 0..m {
        for i in 0..m {
            out[h + m][i] = tau[h][i].clone();
        }
    }

    out
}
