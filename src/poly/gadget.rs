use super::{matrix::PolyMatrix, params::PolyParams, PolyElemOps, PolyOps};

pub trait PolyGadget {
    type Error: std::error::Error + Send + Sync + 'static;
    type M: PolyMatrix;
    type Params;
    fn gadget_vector(params: Self::Params) -> Self::M;
    fn gadget_matrix(params: Self::Params, size: usize) -> Self::M;
    fn decompose(&self) -> Result<Self::M, Self::Error>;
}
