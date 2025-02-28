use super::utils::*;
use super::{Obfuscation, ObfuscationError};
use crate::bgg::*;
use crate::bgg::{eval::*, sampler::*, *};
use crate::poly::gadget::PolyGadgetOps;
use crate::poly::{matrix::*, sampler::*, *};
use crate::utils::*;
use num_traits::{One, Zero};
use rand::Rng;
use rand::RngCore;
use std::sync::Arc;

pub fn obfuscate<T, P, M, G, S, R>(
    poly_op: P,
    matrix_op: M,
    gadget_op: G,
    mut sampler: S,
    rng: &mut R,
    input_size: usize,
) -> Result<Obfuscation<T, P, M>, ObfuscationError>
where
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    G: PolyGadgetOps<T, P, M>,
    S: PolyUniformSampler<T, P, M> + PolyHashSampler<T, P, M> + PolyTrapdoorSampler<T, P, M>,
    R: RngCore,
{
    let hash_key = rng.gen::<[u8; 32]>().to_vec();
    sampler.set_key(&hash_key);
    let sampler = Arc::new(sampler);
    let poly_op = Arc::new(poly_op);
    let matrix_op = Arc::new(matrix_op);
    let gadget_op = Arc::new(gadget_op);
    let deg = poly_op.degree();
    let packed_input_size = ceil_div(input_size, deg);
    let public_data = PublicSampledData::sample(
        matrix_op.clone(),
        gadget_op.clone(),
        sampler.clone(),
        packed_input_size,
    )?;
    let s_bar = sampler
        .sample_uniform::<_, BitDist>(rng, 1, 1)
        .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
    let bgg_encode_sampler = BGGEncodingSampler::new(
        matrix_op.entry(&s_bar, 1, 1).map_err(|e| ObfuscationError::MatrixError(e.to_string()))?,
        sampler.clone(),
        poly_op.clone(),
        matrix_op.clone(),
        gadget_op.clone(),
    );
    let s_init = bgg_encode_sampler.secret_vec.clone();
    let t_bar: PolyMatrix<T, P, M> = sampler
        .sample_uniform::<_, BitDist>(rng, 1, 1)
        .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
    let minus_one_poly = matrix_op.from_poly_vec(vec![poly_op.minus_one()]);
    let t = matrix_op
        .concat_columns(&[t_bar.clone(), minus_one_poly.clone()])
        .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
    let log_q = poly_op.modulus_bits();
    // let a_fhe_bar = sampler
    //     .sample_hash::<_, FinRingDist>(TAG_A_FHE_BAR, 2, 2 * log_q)
    //     .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
    let a_fhe_bar = &public_data.a_fhe_bar;
    let b_fhe = {
        let muled = matrix_op
            .mul(&t_bar, &a_fhe_bar)
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let error = sampler
            .sample_uniform::<_, GaussianDist>(rng, 1, 2 * log_q)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        matrix_op.add(&muled, &error).map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
    };
    // let a_fhe = matrix_op.concat_rows(&[a_fhe_bar, b_fhe]);
    // let fhe_enc =

    // let bgg_pubkeys_input =
    //     bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_INPUT, packed_input_size + 1)?;
    // let bgg_pubkeys_fhe_key = bgg_pubkey_sampler.sample(TAG_BGG_PUBKEY_FHEKEY, 2)?;
    let zero_plaintexts: Vec<Poly<T, P>> = (0..packed_input_size).map(|_| poly_op.zero()).collect();
    let encode_input = bgg_encode_sampler.sample(
        rng,
        &public_data.pubkeys_input,
        &vec![vec![poly_op.one()], zero_plaintexts].concat(),
        true,
    )?;
    let encode_fhe_key = bgg_encode_sampler.sample(
        rng,
        &public_data.pubkeys_fhe_key,
        &matrix_op.to_poly_vec(&t),
        false,
    )?;
    let mut bs = vec![];
    let mut b_trapdoors = vec![];
    for _ in 0..=input_size {
        let (b_0, b_0_trapdoor) =
            sampler.trapdoor(rng).map_err(|e: <S as PolyTrapdoorSampler<T, P, M>>::Error| {
                ObfuscationError::SampleError(e.to_string())
            })?;
        let (b_1, b_1_trapdoor) =
            sampler.trapdoor(rng).map_err(|e: <S as PolyTrapdoorSampler<T, P, M>>::Error| {
                ObfuscationError::SampleError(e.to_string())
            })?;
        let (b_star, b_star_trapdoor) =
            sampler.trapdoor(rng).map_err(|e: <S as PolyTrapdoorSampler<T, P, M>>::Error| {
                ObfuscationError::SampleError(e.to_string())
            })?;
        bs.push((b_0, b_1, b_star));
        b_trapdoors.push((b_0_trapdoor, b_1_trapdoor, b_star_trapdoor));
    }
    let m_b = 2 + log_q;
    let p_init = {
        let s_connect = matrix_op
            .concat_columns(&[s_init.clone(), s_init.clone()])
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let s_b = matrix_op
            .mul(&s_connect, &bs[0].2)
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let error = sampler
            .sample_uniform::<_, GaussianDist>(rng, 1, m_b)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        matrix_op.add(&s_b, &error).map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
    };
    let identity_2 = matrix_op.identity(2, None);
    let u_0 = matrix_op
        .concat_diag(&[identity_2.clone(), public_data.r_0.clone()])
        .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
    let u_1 = matrix_op
        .concat_diag(&[identity_2.clone(), public_data.r_1.clone()])
        .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
    let u_star = {
        let zeros = matrix_op.zero(2, 4);
        let identitys = matrix_op
            .concat_columns(&[identity_2.clone(), identity_2.clone()])
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        matrix_op
            .concat_rows(&[zeros, identitys])
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
    };
    let gadget_2 = gadget_op.gadget_matrix(2);
    let (mut m_preimages, mut n_preimages, mut k_preimages) = (vec![], vec![], vec![]);
    for idx in 0..input_size {
        let (b_next_0, b_next_1, b_next_star) = &bs[idx + 1];
        let (_, _, b_cur_star_trapdoor) = &b_trapdoors[idx];
        let (b_next_0_trapdoor, b_next_1_trapdoor, _) = &b_trapdoors[idx + 1];
        let m_0 = {
            let ub = matrix_op
                .mul(&u_0, &b_next_0)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            sampler
                .preimage(rng, b_cur_star_trapdoor, &ub)
                .map_err(|e| ObfuscationError::SampleError(e.to_string()))?
        };
        let m_1 = {
            let ub: <M as PolyMatrixOps<T, P>>::Matrix = matrix_op
                .mul(&u_1, &b_next_1)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            sampler
                .preimage(rng, b_cur_star_trapdoor, &ub)
                .map_err(|e| ObfuscationError::SampleError(e.to_string()))?
        };
        m_preimages.push((m_0, m_1));

        let ub_star = matrix_op
            .mul(&u_star, &b_next_star)
            .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
        let n_0 = sampler
            .preimage(rng, b_next_0_trapdoor, &ub_star)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        let n_1 = sampler
            .preimage(rng, b_next_1_trapdoor, &ub_star)
            .map_err(|e| ObfuscationError::SampleError(e.to_string()))?;
        n_preimages.push((n_0, n_1));

        let mut ks = vec![];
        for bit in 0..1 {
            let (t_input, t_fhe_key) = if bit == 0 { &public_data.t_0 } else { &public_data.t_1 };
            let at_input = matrix_op
                .mul(&public_data.pubkeys_input[idx].matrix, t_input)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let at_fhe_key = matrix_op
                .mul(&public_data.pubkeys_fhe_key[bit].matrix, t_fhe_key)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let former = matrix_op
                .concat_columns(&[at_input, at_fhe_key])
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let inserted_poly_index = idx / deg;
            let inserted_coeff_index = idx % deg;
            let zero_coeff = <PElem<T> as Zero>::zero();
            let mut coeffs = vec![zero_coeff; deg];
            coeffs[inserted_coeff_index] = <PElem<T> as One>::one();
            let inserted_poly =
                poly_op.from_coeffs(coeffs).map_err(|e: <P as PolyOps<T>>::Error| {
                    ObfuscationError::PolyError(e.to_string())
                })?;
            let inserted_poly_gadget = {
                let zero = poly_op.zero();
                let mut polys = vec![];
                for _ in 0..inserted_poly_index {
                    polys.push(zero.clone());
                }
                polys.push(inserted_poly);
                for _ in inserted_poly_index + 1..packed_input_size {
                    polys.push(zero.clone());
                }
                matrix_op
                    .mul(&matrix_op.from_poly_vec(polys), &gadget_2)
                    .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?
            };
            let a_input_next = matrix_op
                .sub(&public_data.pubkeys_input[idx + 1].matrix, &inserted_poly_gadget)
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let latter = matrix_op
                .concat_columns(&[
                    a_input_next,
                    public_data.pubkeys_fhe_key[idx + 1].matrix.clone(),
                ])
                .map_err(|e| ObfuscationError::MatrixError(e.to_string()))?;
            let k_target = matrix_op.concat_rows(&[former, latter]).map_err(
                |e: <M as PolyMatrixOps<T, P>>::Error| ObfuscationError::MatrixError(e.to_string()),
            )?;
            let trapdoor = if bit == 0 { b_next_0_trapdoor } else { b_next_1_trapdoor };
            let k = sampler.preimage(rng, trapdoor, &k_target).map_err(
                |e: <S as PolyTrapdoorSampler<T, P, M>>::Error| {
                    ObfuscationError::SampleError(e.to_string())
                },
            )?;
            ks.push(k);
        }
        k_preimages.push((ks[0].clone(), ks[1].clone()));
    }
    let obf = Obfuscation {
        hash_key,
        b_fhe,
        encode_input,
        encode_fhe_key,
        p_init,
        m_preimages,
        n_preimages,
        k_preimages,
    };

    Ok(obf)
}
