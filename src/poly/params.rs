pub trait PolyParams: Clone {
    type Modulus: Clone;
    fn modulus(&self) -> Self::Modulus;
    fn modulus_bits(&self) -> usize;
    fn ring_dimension(&self) -> u32;
}
