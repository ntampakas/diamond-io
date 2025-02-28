use super::utils::*;
use super::{Obfuscation, ObfuscationError};
use crate::bgg::*;
use crate::bgg::{eval::*, sampler::*, *};
use crate::poly::gadget::PolyGadgetOps;
use crate::poly::{matrix::*, sampler::*, *};
use crate::utils::*;
use itertools::Itertools;
use num_traits::{One, Zero};
use rand::Rng;
use rand::RngCore;
use std::sync::Arc;

pub fn eval_obf<T, P, M, G, S>(
    poly_op: P,
    matrix_op: M,
    gadget_op: G,
    mut sampler: S,
    obfuscation: Obfuscation<T, P, M>,
    input: &[bool],
) -> Result<Vec<bool>, ObfuscationError>
where
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    G: PolyGadgetOps<T, P, M>,
    S: PolyHashSampler<T, P, M>,
{
    sampler.set_key(&obfuscation.hash_key);
    let sampler = Arc::new(sampler);
    let poly_op = Arc::new(poly_op);
    let matrix_op = Arc::new(matrix_op);
    let gadget_op = Arc::new(gadget_op);
    let deg = poly_op.degree();
    let input_size = input.len();
    let packed_input_size = ceil_div(input_size, deg);
    let public_data = PublicSampledData::sample(
        matrix_op.clone(),
        gadget_op.clone(),
        sampler.clone(),
        packed_input_size,
    )?;
    let (mut ps, mut cs_input, mut cs_fhe_key) = (vec![], vec![], vec![]);
    ps.push(obfuscation.p_init.clone());
    cs_input.push(
        matrix_op
            .concat_columns(
                &obfuscation.encode_input.iter().map(|pubkey| pubkey.vector.clone()).collect_vec(),
            )
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?,
    );
    cs_fhe_key.push(
        matrix_op
            .concat_columns(
                &obfuscation
                    .encode_fhe_key
                    .iter()
                    .map(|pubkey| pubkey.vector.clone())
                    .collect_vec(),
            )
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?,
    );
    let log_q = poly_op.modulus_bits();
    for (idx, input) in input.iter().enumerate() {
        let m =
            if *input { &obfuscation.m_preimages[idx].1 } else { &obfuscation.m_preimages[idx].0 };
        let q = matrix_op
            .mul(&ps[idx], &m)
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let n =
            if *input { &obfuscation.n_preimages[idx].1 } else { &obfuscation.n_preimages[idx].0 };
        let p = matrix_op.mul(&q, &n).map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let k =
            if *input { &obfuscation.k_preimages[idx].1 } else { &obfuscation.k_preimages[idx].0 };
        let v = matrix_op.mul(&q, &k).map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let v_input = matrix_op
            .slice_columns(&v, 0, 2 * log_q * (packed_input_size + 1))
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let v_fhe_key = matrix_op
            .slice_columns(
                &v,
                2 * log_q * (packed_input_size + 1),
                2 * log_q * (packed_input_size + 3),
            )
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let c_input = {
            // let encode =
            let t = if *input { &public_data.t_1.0 } else { &public_data.t_0.0 };
            let muled = matrix_op
                .mul(&cs_input[idx], &t)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            matrix_op
                .add(&muled, &v_input)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
        };
        let c_fhe_key = {
            let t = if *input { &public_data.t_1.1 } else { &public_data.t_0.1 };
            let muled = matrix_op
                .mul(&cs_fhe_key[idx], &t)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            matrix_op
                .add(&muled, &v_fhe_key)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
        };
        ps.push(p);
        cs_input.push(c_input);
        cs_fhe_key.push(c_fhe_key);
    }
    Ok(vec![])
}
