use num_bigint::BigUint;

pub trait PolyParams: Clone {
    fn modulus(&self) -> BigUint;
    fn modulus_bits(&self) -> usize {
        self.modulus().bits() as usize
    }
    fn ring_dimension(&self) -> u32;
}
