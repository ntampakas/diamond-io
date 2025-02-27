use std::marker::PhantomData;

use digest::OutputSizeUser;
use num_bigint::BigUint;

use crate::poly::{
    dcrt::{DCRTPoly, DCRTPolyMatrix, FieldElement},
    PolyParams, Polynomial, PolynomialMatrix,
};

use super::PolyHashSamplerTrait;

pub enum PolyHashDistType {
    FinRingDist,
    BitDist,
}

pub struct PolyHashSampler<P, M, D>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
    D: digest::Digest,
{
    key: [u8; 32],
    dist_type: PolyHashDistType,
    params: P::Params,
    _phantom_p: PhantomData<P>,
    _phantom_m: PhantomData<M>,
    _phantom_d: PhantomData<D>,
}

impl<D> PolyHashSampler<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, D>
where
    D: digest::Digest,
{
    pub fn new(key: [u8; 32], dist_type: PolyHashDistType, params: PolyParams) -> Self {
        Self {
            key,
            dist_type,
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
        let hash_output_size = <D as digest::Digest>::output_size() * 8;
        let n = self.params.get_ring_dimension() as usize; // TODO: fix type
        let q = self.params.get_modulus() as usize; // TODO: fix type
        let field_elements = match self.dist_type {
            PolyHashDistType::FinRingDist => {
                // * index = number of hashes to be performed = ceil( (nrow * ncol * n * ilog2(q)) / (hash_output_size) )
                // * bits = number of resulting bits from hashing ops = hash_output_size * index
                // * field_elements = number of field elements sampled = (bits / ilog2(q)) which is always greater than or equal to nrow * ncol * n
                let log2q = q.ilog2() as usize;
                let index = (nrow * ncol * n * log2q).div_ceil(hash_output_size);
                let mut bits = Vec::with_capacity(hash_output_size * index);
                let mut field_elements = Vec::with_capacity((index * hash_output_size) / log2q);
                for i in 0..(index) {
                    //  H ( key || tag || i )
                    let mut hasher = D::new();
                    let mut combined = Vec::with_capacity(self.key.len() + tag.as_ref().len() + 1);
                    combined.extend_from_slice(&self.key);
                    combined.extend_from_slice(tag.as_ref());
                    combined.push(i as u8);
                    hasher.update(&combined);
                    for &byte in hasher.finalize().iter() {
                        for bit_index in 0..8 {
                            let bit = (byte >> bit_index) & 1;
                            bits.push(bit);
                        }
                    }
                }
                // From bits to field elements
                let mut offset = 0;
                for _ in 0..(bits.len() / log2q) {
                    let value_bits = &bits[offset..offset + log2q];
                    let value = BigUint::from_radix_be(value_bits, 2).unwrap();
                    offset += log2q;
                    let fe = FieldElement::new(value, BigUint::from(q));
                    field_elements.push(fe);
                }
                field_elements
            }
            PolyHashDistType::BitDist => {
                // * index = number of hashes to be performed = ceil( (nrow * ncol * n) / (hash_output_size) )
                // * bits = number of resulting bits from hashing ops = hash_output_size * index
                // * field_elements = number of field elements sampled = bits
                // which is always greater than or equal to nrow * ncol * n
                let index = (nrow * ncol * n).div_ceil(hash_output_size);
                let mut field_elements = Vec::with_capacity(hash_output_size * index);
                for i in 0..index {
                    //  H ( key || tag || i )
                    let mut hasher = D::new();
                    let mut combined = Vec::with_capacity(self.key.len() + tag.as_ref().len() + 1);
                    combined.extend_from_slice(&self.key);
                    combined.extend_from_slice(tag.as_ref());
                    combined.push(i as u8);
                    hasher.update(&combined);
                    for &byte in hasher.finalize().iter() {
                        for bit_index in 0..8 {
                            let bit = (byte >> bit_index) & 1;
                            field_elements.push(FieldElement::new(bit as u64, q as u64));
                        }
                    }
                }
                field_elements
            }
        };

        // Check if we have enough field elements (not sure if this will ever happen)
        if field_elements.len() < nrow * ncol * n {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                "Not enough field elements to sample hash",
            ));
        }

        // From field elements to nrow * ncol polynomials
        let total_poly = nrow * ncol;
        let mut offset = 0;
        let mut all_polys = Vec::with_capacity(total_poly);
        for _ in 0..total_poly {
            let coeffs = &field_elements[offset..offset + n];
            offset += n;
            let poly = DCRTPoly::from_coeffs(&self.params, coeffs)?;
            all_polys.push(poly);
        }

        // From polynomials to matrix such that the first row of the matrixcontains the first ncol polynomials
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
        let mut sampler = PolyHashSampler::<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, Keccak256>::new(
            key,
            PolyHashDistType::BitDist,
            params,
        );
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(tag, nrow, ncol);

        let new_key = [1u8; 32];
        sampler.set_key(new_key);
        sampler.expose_key();

        assert!(matrix_result.is_ok(), "sample_hash returned an error: {:?}", matrix_result.err());
        let matrix = matrix_result.unwrap();
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }

    #[test]
    fn test_poly_hash_sampler_fin_ring_dist() {
        let key = [0u8; 32];
        let params = PolyParams::new(16, 4, 51);
        let mut sampler = PolyHashSampler::<DCRTPoly, DCRTPolyMatrix<DCRTPoly>, Keccak256>::new(
            key,
            PolyHashDistType::FinRingDist,
            params,
        );
        let nrow = 100;
        let ncol = 300;
        let tag = b"MyTag";
        let matrix_result = sampler.sample_hash(tag, nrow, ncol);

        let new_key = [1u8; 32];
        sampler.set_key(new_key);
        sampler.expose_key();

        assert!(matrix_result.is_ok(), "sample_hash returned an error: {:?}", matrix_result.err());
        let matrix = matrix_result.unwrap();
        assert_eq!(matrix.row_size(), nrow, "Matrix row count mismatch");
        assert_eq!(matrix.col_size(), ncol, "Matrix column count mismatch");
    }
}
