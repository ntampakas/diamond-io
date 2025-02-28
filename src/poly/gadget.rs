use super::matrix::PolyMatrix;

pub trait PolyGadget {
    type Error: std::error::Error + Send + Sync + 'static;
    type M: PolyMatrix;
    type Params;
    // TODO: do we need on interface level?
    fn gadget_vector(params: &Self::Params) -> Self::M;
    fn gadget_matrix(params: &Self::Params, size: usize) -> Self::M;
    fn decompose(&self) -> Result<Self::M, Self::Error>;
}
