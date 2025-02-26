use std::marker::PhantomData;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix},
    Polynomial, PolynomialMatrix,
};

use super::PolyHashSamplerTrait;

pub struct PolyHashSampler<P, M>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
{
    _phantom_p: PhantomData<P>,
    _phantom_m: PhantomData<M>,
}

impl PolyHashSamplerTrait<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
    for PolyHashSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>>
{
    type Error = std::io::Error;
    type Key = [u8; 32];

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        _tag: B,
        _nrow: usize,
        _ncol: usize,
    ) -> Result<DCRTPolyMatrix<DCRTPoly>, Self::Error> {
        todo!()
    }

    fn set_key(&mut self, _key: Self::Key) {
        todo!()
    }

    fn expose_key(&self) -> Self::Key {
        todo!()
    }
}
