pub mod gadget;
pub mod matrix;
pub mod sampler;
use num_traits::{One, Zero};
use phantom_zone_math::modulus::{Elem, ElemFrom, ElemOps, ElemTo};
use phantom_zone_math::util::serde::Serde;
use std::fmt::Debug;
use std::hash::Hash;

pub trait PolyElemModulus {
    fn modulus_bits(&self) -> usize;
}

pub trait PolyDegree<T: PolyElemOps>: PolyElemModulus {
    fn degree(&self) -> usize;
}

pub type PElem<T> = <T as PolyElemOps>::PElem;

pub trait PolyElemOps: Clone + Debug + Send + Sync + PolyElemModulus {
    type Error: std::error::Error + Send + Sync + 'static;
    type PElem: 'static
        + Copy
        + Debug
        + Default
        + Eq
        + Ord
        + Hash
        + Send
        + Sync
        + Serde
        + Zero
        + One;
}

pub type Poly<T, P> = <P as PolyOps<T>>::Poly;

pub trait PolyOps<T: PolyElemOps>: Clone + Debug + Send + Sync + PolyDegree<T> {
    type Error: std::error::Error + Send + Sync + 'static;
    type Poly: Debug + Clone;
    fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    fn from_coeffs(&self, coeffs: Vec<PElem<T>>) -> Result<Self::Poly, Self::Error>;
    fn from_const(&self, constant: &T) -> Result<Self::Poly, Self::Error>;
    fn zero(&self) -> Self::Poly;
    fn one(&self) -> Self::Poly;
    fn minus_one(&self) -> Self::Poly;
    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn neg(&self, a: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let minus_b = self.neg(b)?;
        self.add(a, &minus_b)
    }
    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error>;
}
