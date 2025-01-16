mod operations;
mod parameters;
mod utils;

pub use operations::{poly_add, poly_sub, pub_key_gen, Ciphertext};
pub use parameters::Parameters;
pub use utils::{print_matrix_ring, print_vector, print_vector_ring};

pub struct BggRlwe {
    params: Parameters,
    public_key: Vec<Vec<Vec<u64>>>,
}

impl BggRlwe {
    pub fn new(log_ring_size: usize, k: usize, ell: usize) -> Self {
        let params = Parameters::new(log_ring_size, k, ell);
        let public_key = operations::pub_key_gen(&params);
        Self { params, public_key }
    }

    /// Encode an attribute vector into a ciphertext
    ///
    /// # Returns
    /// - `ct_inner`: The inner component ciphertext `(1, x) ⊗ G + B`
    /// - `ct_in`: The entire ciphertext `((1, x) ⊗ G + B) s + e_0`
    pub fn encode(&self, x: Vec<u64>) -> Ciphertext {
        Ciphertext::new(&self.public_key, &self.params, &x)
    }

    // pub fn m_eval_c(&self, C) TO ADD - takes a circuit in Bristol format as input and returns the pub key evaluated on C
    // pub fn m_eval_c_x(&self, C) TO ADD - takes a circuit in Bristol format as input and returns the ciphertext evaluated on C
}
