pub mod dcrt_matrix;
pub mod dcrt_poly;
pub mod field_element;
// pub mod gadget;
pub mod matrix;
pub mod params;
pub mod sampler;
use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Mul, MulAssign, Neg},
};

use num_traits::{One, Zero};
use params::Params;

// use phantom_zone_math::modulus::Elem;
// pub type PElem<T> = Elem<T>;

// pub trait PolyElemModulus {
//     fn modulus_bits(&self) -> usize;
// }

pub trait PElem:
    'static
    + Copy
    + Debug
    + Default
    + Eq
    + Ord
    // + Hash
    + Send
    + Sync
    // + Serialize
    // + Deserialize<'static>
    // + Zero TODO: do we need this?
    // + One TODO: do we need this?
    + Add
    + Mul
{
}

// pub trait PolyElemOps:
//     ElemOps + ElemFrom<u64> + ElemFrom<u32> + ElemFrom<u8> + ElemFrom<bool> + PolyElemModulus
// {
//     type Error: std::error::Error + Send + Sync + 'static;
// }

// pub trait PolyBitOps: PolyElemOps + ElemTo<u64> + ElemTo<u32> + ElemTo<u8> + ElemTo<bool> {
//     fn modulus_bits(&self) -> usize {
//         1
//     }
// }

// pub trait PolyGaussOps: PolyElemOps {
//     fn gaussian_param(&self) -> f64;
// }

// pub type Poly<T, P> = <P as PolyOps<T>>::Poly;

// pub trait PolyDegree {
//     fn degree(&self) -> usize;
// }

/// Describes the common interface polynomials
pub trait Polynomial:
    Sized + Clone + Debug + PartialEq + Eq + Add + AddAssign + Mul + MulAssign + Neg + Zero + One
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Elem: PElem;
    // fn degree(&self) -> usize;
    // fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    // fn from_coeffs(coeffs: &[PElem<T>]) -> Result<Self::Poly, Self::Error>;
    fn from_const(params: &Params, constant: &Self::Elem) -> Result<Self, Self::Error>; // TODO: replace with u64 with T
    fn const_zero(params: &Params) -> Self;
    fn const_one(params: &Params) -> Self;
    // fn const_minus_one(params: &Params) -> Result<Self, Self::Error>;
    // fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error>;
}
