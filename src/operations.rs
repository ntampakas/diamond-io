use crate::{
    parameters::Parameters,
    utils::{empty_matrix_ring, empty_ring_element, empty_vector_ring},
};
use phantom_zone_math::{
    prelude::{ElemFrom, ModulusOps},
    ring::{PrimeRing, RingOps},
};

pub fn bit_decompose(params: &Parameters, bu: &Vec<Vec<u64>>) -> Vec<Vec<Vec<u64>>> {
    let ring = params.ring();
    let m = *params.m();
    let ring_size = ring.ring_size();
    // Create a matrix of dimension m × m, where each element is a binary polynomial
    let mut tau = empty_matrix_ring(ring, m, m);

    // For each row h in the output matrix
    for h in 0..m {
        // For each column i in the output matrix
        for i in 0..m {
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

pub fn poly_add(ring: &PrimeRing, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    let mut c = empty_ring_element(ring);
    for i in 0..a.len() {
        let elem = ring.add(&a[i], &b[i]);
        c[i] = elem;
    }
    c
}

pub fn poly_sub(ring: &PrimeRing, a: &Vec<u64>, b: &Vec<u64>) -> Vec<u64> {
    assert_eq!(a.len(), b.len());
    let mut c = empty_ring_element(ring);
    for i in 0..a.len() {
        let elem = ring.sub(&a[i], &b[i]);
        c[i] = elem;
    }
    c
}

pub fn gen_identity_matrix_to_scalar(
    ring: &PrimeRing,
    m: usize,
    scalar: u64,
) -> Vec<Vec<Vec<u64>>> {
    let mut identity_matrix = empty_matrix_ring(ring, m, m);

    for i in 0..m {
        identity_matrix[i][i][0] = ring.elem_from(scalar);
    }

    identity_matrix
}

pub fn vec_vec_dot_product(
    ring: &PrimeRing,
    vec_a: &Vec<Vec<u64>>,
    vec_b: &Vec<Vec<u64>>,
) -> Vec<u64> {
    assert_eq!(vec_a.len(), vec_b.len(),);
    let mut out = empty_ring_element(ring);
    for i in 0..vec_a.len() {
        let mut scratch = ring.allocate_scratch(1, 2, 0);
        let mut scratch = scratch.borrow_mut();
        let product = ring.take_poly(&mut scratch);
        ring.poly_mul(product, &vec_a[i], &vec_b[i], scratch.reborrow());
        out = poly_add(ring, &out, &product.to_vec());
    }
    out
}

pub fn vec_mat_mul(
    ring: &PrimeRing,
    vec: &Vec<Vec<u64>>,
    mat: &Vec<Vec<Vec<u64>>>,
) -> Vec<Vec<u64>> {
    let len = vec.len();
    let mat_rows = mat.len();
    let mat_cols = mat[0].len();
    assert_eq!(len, mat_rows);
    let mut out = empty_vector_ring(ring, mat_cols);

    for i in 0..mat_cols {
        let col_i = mat.iter().map(|row| row[i].clone()).collect::<Vec<_>>();
        out[i] = vec_vec_dot_product(ring, vec, &col_i);
    }

    out
}

pub fn mat_mat_add(
    ring: &PrimeRing,
    mat_a: &Vec<Vec<Vec<u64>>>,
    mat_b: &Vec<Vec<Vec<u64>>>,
) -> Vec<Vec<Vec<u64>>> {
    assert_eq!(mat_a.len(), mat_b.len());
    assert_eq!(mat_a[0].len(), mat_b[0].len());
    assert_eq!(mat_a.len(), mat_b.len());
    assert_eq!(mat_a[0].len(), mat_b[0].len());

    let mat_a_rows = mat_a.len();
    let mat_a_cols = mat_a[0].len();

    let mut out = empty_matrix_ring(ring, mat_a_rows, mat_a_cols);

    for i in 0..mat_a.len() {
        for j in 0..mat_a[0].len() {
            out[i][j] = poly_add(ring, &mat_a[i][j], &mat_b[i][j]);
        }
    }
    out
}

pub fn vec_vec_add(
    ring: &PrimeRing,
    vec_a: &Vec<Vec<u64>>,
    vec_b: &Vec<Vec<u64>>,
) -> Vec<Vec<u64>> {
    assert_eq!(vec_a.len(), vec_b.len());
    let mut out = empty_vector_ring(ring, vec_a.len());
    for i in 0..vec_a.len() {
        out[i] = poly_add(ring, &vec_a[i], &vec_b[i]);
    }
    out
}
/// Vertically concatenate two matrices of ring elements of equal size
pub fn mat_vert_concat(
    ring: &PrimeRing,
    mat_top: &Vec<Vec<Vec<u64>>>,
    mat_bottom: &Vec<Vec<Vec<u64>>>,
) -> Vec<Vec<Vec<u64>>> {
    assert_eq!(mat_top.len(), mat_bottom.len());
    assert_eq!(mat_top[0].len(), mat_bottom[0].len(),);
    assert_eq!(mat_top[0][0].len(), mat_bottom[0][0].len(),);

    let rows = mat_top.len();
    let cols = mat_top[0].len();

    // Create output matrix with combined number of rows
    let mut out = empty_matrix_ring(ring, 2 * rows, cols);

    // Copy top matrix
    for i in 0..rows {
        for j in 0..cols {
            out[i][j] = mat_top[i][j].clone();
        }
    }

    // Copy bottom matrix
    for i in 0..rows {
        for j in 0..cols {
            out[i + rows][j] = mat_bottom[i][j].clone();
        }
    }

    out
}

/// Horizontally concatenate two matrices of ring elements of equal size
pub fn mat_horiz_concat(
    ring: &PrimeRing,
    mat_left: &Vec<Vec<Vec<u64>>>,
    mat_right: &Vec<Vec<Vec<u64>>>,
) -> Vec<Vec<Vec<u64>>> {
    assert_eq!(mat_left.len(), mat_right.len());
    assert_eq!(mat_left[0].len(), mat_right[0].len());
    assert_eq!(mat_left[0][0].len(), mat_right[0][0].len());

    let rows = mat_left.len();
    let cols = mat_right[0].len();

    // Create output matrix with combined number of columns
    let mut out = empty_matrix_ring(ring, rows, 2 * cols);

    // Copy left matrix
    for i in 0..rows {
        for j in 0..cols {
            out[i][j] = mat_left[i][j].clone();
        }
    }

    // Copy right matrix
    for i in 0..rows {
        for j in 0..cols {
            out[i][j + cols] = mat_right[i][j].clone();
        }
    }

    out
}

pub fn vec_horiz_concat(vec_left: &Vec<Vec<u64>>, vec_right: &Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    assert_eq!(vec_left.len(), vec_right.len());
    let mut out = vec_left.clone();
    out.extend(vec_right.clone());
    out
}

#[cfg(test)]
mod tests {
    use phantom_zone_math::{prelude::ModulusOps, ring::RingOps};

    use crate::{
        operations::{bit_decompose, poly_add},
        parameters::Parameters,
        pub_key::PublicKey,
    };
    #[test]
    fn test_bit_decompose() {
        let params = Parameters::new(12, 51, 4);
        let pub_key = PublicKey::new(params);
        let b1 = &pub_key.b()[1];
        let ring = pub_key.params().ring();
        let m = *pub_key.params().m();
        let g = pub_key.params().g();
        let tau = bit_decompose(pub_key.params(), b1);

        // Reconstruct the original input by multiplying tau with G
        let mut reconstructed = vec![vec![ring.zero(); ring.ring_size()]; m];

        // For each column i of the output
        for i in 0..m {
            // For each row h of tau
            for h in 0..m {
                // Multiply tau[h][i] by g[h] and add to the result
                let mut scratch = ring.allocate_scratch(1, 2, 0);
                let mut scratch = scratch.borrow_mut();
                let product = ring.take_poly(&mut scratch);
                ring.poly_mul(product, &tau[h][i], &g[h], scratch.reborrow());
                reconstructed[i] = poly_add(ring, &reconstructed[i], &product.to_vec());
            }
        }

        for i in 0..m {
            assert_eq!(b1[i], reconstructed[i]);
        }
    }
}
