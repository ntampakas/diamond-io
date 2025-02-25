pub mod eval;
pub mod obf;

use crate::bgg::{BggEncoding, BggError, BggPublicKey};
use crate::poly::gadget::PolyGadgetOps;
use crate::poly::{matrix::*, sampler::*, *};

#[derive(Debug, Clone)]
pub struct Obfuscation<T, P, M, S>
where
    T: PolyElemOps,
    P: PolyOps<T>,
    M: PolyMatrixOps<T, P>,
    S: PolyUniformSampler<T, P, M, BitDist>
        + PolyUniformSampler<T, P, M, GaussianDist>
        + PolyUniformSampler<T, P, M, FinRingDist>
        + PolyUniformSampler<T, P, M, BitDist>,
{
    pub sampler: S,
    pub fhe_enc: PolyMatrix<T, P, M>,
    pub input_encode: BggEncoding<T, P, M>,
    pub fhe_key_encode: BggPublicKey<T, P, M>,
    pub init_p: PolyMatrix<T, P, M>,
    pub m_preimages: PolyMatrix<T, P, M>,
    pub n_preimages: PolyMatrix<T, P, M>,
    pub k_preimages: PolyMatrix<T, P, M>,
    pub final_preimage: PolyMatrix<T, P, M>,
}
