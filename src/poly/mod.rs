pub mod dcrt_poly;
// pub mod gadget;
// pub mod matrix;
pub mod params;
pub mod sampler;
use params::Params;
use phantom_zone_math::modulus::Elem;
pub type PElem<T> = Elem<T>;

pub trait PolyElemModulus {
    fn modulus_bits(&self) -> usize;
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

pub trait PolyOps {
    type Error: std::error::Error + Send + Sync + 'static; // TODO: we should print UniquePtr
    type Poly;
    // fn degree(&self) -> usize;
    // fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    // fn from_coeffs(coeffs: &[PElem<T>]) -> Result<Self::Poly, Self::Error>;
    fn from_const(params: &Params, constant: &u64) -> Result<Self::Poly, Self::Error>; // TODO: replace with u64 with T
    // fn zero(&self) -> Result<Self::Poly, Self::Error>;
    // fn one(&self) -> Result<Self::Poly, Self::Error>;
    // fn minus_one(&self) -> Result<Self::Poly, Self::Error>;
    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    // fn neg(&self, a: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    // fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
    //     let minus_b = self.neg(b)?;
    //     self.add(a, &minus_b)
    // }
    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    // fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error>;
}
