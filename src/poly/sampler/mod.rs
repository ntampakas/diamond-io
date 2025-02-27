use super::{Poly, PolynomialMatrix};

pub mod hash;
pub mod trapdoor;
pub mod uniform;

pub trait PolyHashSampler<P, M, D>
where
    P: Poly,
    M: PolynomialMatrix<P>,
    D: digest::Digest,
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Key; // TODO: what to do with this type?

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        rows: usize,
        columns: usize,
    ) -> Result<M, Self::Error>;

    fn set_key(&mut self, key: Self::Key);

    fn expose_key(&self) -> Self::Key;
}

pub trait UniformSampler<P, M>
where
    P: Poly,
    M: PolynomialMatrix<P>,
{
    type Error: std::error::Error + Send + Sync + 'static;

    fn sample_uniform(&self, rows: usize, columns: usize) -> Result<M, Self::Error>;
}
