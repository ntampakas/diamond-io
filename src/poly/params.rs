use num_bigint::BigUint;

pub trait PolyParams: Clone {
    /// Returns the modulus value `q` used for polynomial coefficients in the ring `Z_q[x]/(x^n - 1)`.
    fn modulus(&self) -> BigUint;
    /// Fewest bits necessary to represent the modulus value `q`.
    fn modulus_bits(&self) -> usize {
        self.modulus().bits() as usize
    }
    /// Returns the integer `n` that specifies the size of the polynomial ring used in this
    /// polynomial. Specifically, this is the degree parameter for the ring `Z_q[x]/(x^n - 1)`.
    fn ring_dimension(&self) -> u32;
}
