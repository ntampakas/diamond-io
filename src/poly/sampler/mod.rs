use super::{Polynomial, PolynomialMatrix};

pub mod hash;
pub mod trapdoor;
pub mod uniform;

pub trait PolyHashSamplerTrait<P, M>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
{
    type Error: std::error::Error + Send + Sync + 'static;
    type Key;

    fn sample_hash<B: AsRef<[u8]>>(
        &self,
        tag: B,
        rows: usize,
        columns: usize,
    ) -> Result<M, Self::Error>;

    fn set_key(&mut self, key: Self::Key);

    fn expose_key(&self) -> Self::Key;
}

pub trait MatrixUniformSamplerTrait<P, M>
where
    P: Polynomial,
    M: PolynomialMatrix<P>,
{
    type Error: std::error::Error + Send + Sync + 'static;

    fn sample_uniform(&self, rows: usize, columns: usize) -> Result<M, Self::Error>;
}
