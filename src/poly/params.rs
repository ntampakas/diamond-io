use num_bigint::BigUint;

// pub trait PolyElemParams: Clone {
//     fn modulus(&self) -> BigUint;
// }

pub trait PolyParams: Clone {
    fn modulus(&self) -> BigUint;

    fn ring_dimension(&self) -> u32;
}
