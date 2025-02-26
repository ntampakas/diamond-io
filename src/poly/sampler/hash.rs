use std::marker::PhantomData;

use digest::OutputSizeUser;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix, FieldElement},
    PolyParams, Polynomial, PolynomialMatrix,
};

use super::PolyHashSamplerTrait;

pub struct PolyHashSampler<P, M, D>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
    D: digest::Digest,
{
    key: [u8; 32],
    params: P::Params,
    _phantom_p: PhantomData<P>,
    _phantom_m: PhantomData<M>,
    _phantom_d: PhantomData<D>,
}

impl<D> PolyHashSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, D>
where
    D: digest::Digest,
{
    pub fn new(key: [u8; 32], params: PolyParams) -> Self {
        Self {
            key,
            params,
            _phantom_p: PhantomData,
            _phantom_m: PhantomData,
            _phantom_d: PhantomData,
        }
    }
}

impl<D> PolyHashSamplerTrait<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, D>
    for PolyHashSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, D>
where
    D: OutputSizeUser + digest::Digest,
{
    type Error = std::io::Error;
    type Key = [u8; 32];

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        nrow: usize,
        ncol: usize,
    ) -> Result<DCRTPolyMatrix<DCRTPoly>, Self::Error> {
        // 1) Concatenate [ key || tag || index ]
        let hash_output_size = <D as digest::Digest>::output_size() * 8;
        // todo: get from params
        let n = 10;
        // degree of polynomial
        let q = 10;
        let index = (nrow * ncol * q).div_ceil(hash_output_size);

        let mut bits = Vec::with_capacity(hash_output_size * index);
        for i in 0..(index) {
            let mut hasher = D::new();
            let mut combined = Vec::with_capacity(self.key.len() + tag.as_ref().len() + 1);
            combined.extend_from_slice(&self.key);
            combined.extend_from_slice(tag.as_ref());
            combined.push(i as u8);
            hasher.update(&combined);
            let result = hasher.finalize();

            for &byte in result.iter() {
                for bit_index in 0..8 {
                    let bit = (byte >> bit_index) & 1;
                    bits.push(FieldElement::new(bit as u64, q as u64));
                }
            }
        }

        let total_poly = nrow * ncol;
        let mut offset = 0;
        let mut all_polys = Vec::with_capacity(total_poly);
        for _ in 0..total_poly {
            let coeffs = &bits[offset..offset + n];
            offset += n;
            let poly = DCRTPoly::from_coeffs(&self.params, coeffs)?;
            all_polys.push(poly);
        }
        let mut matrix_inner = Vec::with_capacity(nrow);
        let mut poly_iter = all_polys.into_iter();
        for _ in 0..nrow {
            let row_polys: Vec<DCRTPoly> = poly_iter.by_ref().take(ncol).collect();
            matrix_inner.push(row_polys);
        }

        Ok(DCRTPolyMatrix::from_poly_vec(&self.params, matrix_inner))
    }

    fn set_key(&mut self, key: Self::Key) {
        self.key = key
    }

    fn expose_key(&self) -> Self::Key {
        self.key
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use keccak_asm::Keccak256;

    #[test]
    fn test_poly_hash_sampler() {
        let key = [0u8; 32];
        let params = PolyParams::new(16, 4, 51);
        let mut sampler =
            PolyHashSampler::<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, Keccak256>::new(key, params);
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(tag, nrow, ncol);

        let new_key = [0u8; 32];
        sampler.set_key(new_key);
        sampler.expose_key();

        assert!(matrix_result.is_ok(), "sample_hash returned an error: {:?}", matrix_result.err());
        let matrix = matrix_result.unwrap();
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }
}
