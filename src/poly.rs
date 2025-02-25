pub mod gadget;
pub mod matrix;
pub mod sampler;
use phantom_zone_math::modulus::{Elem, ElemFrom, ElemOps, ElemTo};
use phantom_zone_math::util::serde;
use std::fmt::Debug;

pub type PElem<T> = Elem<T>;

pub trait PolyElemModulus {
    fn modulus_bits(&self) -> usize;
}
pub trait PolyElemOps:
    ElemOps + ElemFrom<u64> + ElemFrom<u32> + ElemFrom<u8> + ElemFrom<bool> + PolyElemModulus
{
    type Error: std::error::Error + Send + Sync + 'static;
}

// pub trait PolyBitOps: PolyElemOps + ElemTo<u64> + ElemTo<u32> + ElemTo<u8> + ElemTo<bool> {
//     fn modulus_bits(&self) -> usize {
//         1
//     }
// }

// pub trait PolyGaussOps: PolyElemOps {
//     fn gaussian_param(&self) -> f64;
// }

pub type Poly<T, P> = <P as PolyOps<T>>::Poly;

pub trait PolyDegree<T: PolyElemOps>: PolyElemModulus {
    fn degree(&self) -> usize;
}

pub trait PolyOps<T: PolyElemOps>: PolyDegree<T> {
    type Error: std::error::Error + Send + Sync + 'static;
    type Poly: Debug + Clone;
    fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    fn from_coeffs(&self, coeffs: &[PElem<T>]) -> Result<Self::Poly, Self::Error>;
    fn from_const(&self, constant: &T) -> Result<Self::Poly, Self::Error>;
    fn zero(&self) -> Result<Self::Poly, Self::Error>;
    fn one(&self) -> Result<Self::Poly, Self::Error>;
    fn minus_one(&self) -> Result<Self::Poly, Self::Error>;
    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn neg(&self, a: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let minus_b = self.neg(b)?;
        self.add(a, &minus_b)
    }
    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error>;
}
