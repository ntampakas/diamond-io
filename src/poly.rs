pub mod matrix;
pub mod sampler;
use phantom_zone_math::modulus::{Elem, ElemFrom, ElemOps, ElemTo};
use phantom_zone_math::util::serde;

pub type PElem<T> = Elem<T>;
pub trait PolyElemOps:
    ElemOps + ElemFrom<u64> + ElemFrom<u32> + ElemFrom<u8> + ElemFrom<bool>
{
    fn modulus_bits(&self) -> usize;
}

pub trait PolyBitOps: PolyElemOps + ElemTo<u64> + ElemTo<u32> + ElemTo<u8> + ElemTo<bool> {
    fn modulus_bits(&self) -> usize {
        1
    }
}

pub trait PolyGaussOps: PolyElemOps {
    fn gaussian_param(&self) -> f64;
}

pub type Poly<T, P> = <P as PolyOps<T>>::Poly;
pub trait PolyOps<T: PolyElemOps> {
    type Error: std::error::Error;
    type Poly;
    fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    fn from_coeffs(coeffs: &[PElem<T>]) -> Result<Self::Poly, Self::Error>;
    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn decompose(&self, poly: &Self::Poly, base: u8) -> Result<Vec<Self::Poly>, Self::Error>;
}

// #[cfg(test)]
// pub mod tests {
//     use super::*;

//     #[derive(Debug, Clone)]
//     pub struct MockElemOps;

//     #[derive(Debug, Clone)]
//     pub struct MockElem(u64);

//     impl ElemOps for MockElemOps {
//         type Elem = MockElem;
//     }

//     impl ElemFrom<u64> for MockElemOps {
//         fn elem_from(&self, value: u64) -> <Self as ElemOps>::Elem {
//             MockElem(value)
//         }
//     }

//     impl ElemFrom<u32> for MockElemOps {
//         fn elem_from(&self, value: u32) -> <Self as ElemOps>::Elem {
//             MockElem(value as u64)
//         }
//     }

//     impl ElemFrom<u8> for MockElemOps {
//         fn elem_from(&self, value: u8) -> <Self as ElemOps>::Elem {
//             MockElem(value as u64)
//         }
//     }

//     impl ElemFrom<bool> for MockElemOps {
//         fn elem_from(&self, value: bool) -> <Self as ElemOps>::Elem {
//             MockElem(value as u64)
//         }
//     }

//     impl ElemTo<u64> for MockElemOps {
//         fn elem_to(&self, elem: <Self as ElemOps>::Elem) -> u64 {
//             elem.0
//         }
//     }

//     impl ElemTo<u32> for MockElemOps {
//         fn elem_to(&self, elem: <Self as ElemOps>::Elem) -> u32 {
//             elem.0 as u32
//         }
//     }

//     impl ElemTo<u8> for MockElemOps {
//         fn elem_to(&self, elem: <Self as ElemOps>::Elem) -> u8 {
//             elem.0 as u8
//         }
//     }

//     impl ElemTo<bool> for MockElemOps {
//         fn elem_to(&self, elem: <Self as ElemOps>::Elem) -> bool {
//             elem.0 != 0
//         }
//     }

//     impl PolyElemOps for MockElemOps {
//         fn modulus_bits(&self) -> usize {
//             64
//         }
//     }

//     impl PolyBitOps for MockElemOps {}

//     impl PolyGaussOps for MockElemOps {
//         fn gaussian_param(&self) -> f64 {
//             3.0
//         }
//     }

//     #[derive(Debug, Clone)]
//     pub struct MockPoly {
//         coeffs: Vec<PElem<MockElemOps>>,
//     }

//     impl Poly<MockElemOps> for MockPoly {
//         fn coeffs(&self) -> &[PElem<MockElemOps>] {
//             &self.coeffs
//         }

//         fn from_coeffs(coeffs: Vec<PElem<MockElemOps>>) -> Self {
//             Self { coeffs }
//         }
//     }
// }
