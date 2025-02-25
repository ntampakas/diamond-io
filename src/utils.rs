use phantom_zone_math::{
    prelude::ModulusOps,
    ring::{PrimeRing, RingOps},
};

pub fn ceil_div(a: usize, b: usize) -> usize {
    a.div_ceil(b)
}

/// Print a ring element
pub fn print_ring_element(label: &str, ring_el: &[u64]) {
    print!("{} [", label);
    for (k, &val) in ring_el.iter().enumerate() {
        if k > 0 {
            print!(", ");
        }
        print!("{}", val);
    }
    println!("]");
}

/// Print a matrix of ring elements
pub fn print_matrix_ring(label: &str, matrix: &[Vec<Vec<u64>>]) {
    println!("\n{}", label,);

    for (i, row) in matrix.iter().enumerate() {
        for (j, col) in row.iter().enumerate() {
            print!("r{}c{}: ", i, j);
            print_ring_element("", col);
        }
    }
}

/// Print a vector of ring elements
pub fn print_vector_ring(label: &str, vec: &[Vec<u64>]) {
    println!("\n{}", label);
    for (i, inner_vec) in vec.iter().enumerate() {
        print!("{}[{}]: ", label, i);
        print_ring_element("", inner_vec);
    }
}

pub fn empty_matrix_ring(ring: &PrimeRing, rows: usize, cols: usize) -> Vec<Vec<Vec<u64>>> {
    vec![vec![vec![ring.zero(); ring.ring_size()]; cols]; rows]
}

pub fn empty_vector_ring(ring: &PrimeRing, cols: usize) -> Vec<Vec<u64>> {
    vec![vec![ring.zero(); ring.ring_size()]; cols]
}

pub fn empty_ring_element(ring: &PrimeRing) -> Vec<u64> {
    vec![ring.zero(); ring.ring_size()]
}
