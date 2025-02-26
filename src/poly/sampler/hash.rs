// pub trait PolyHashSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D: DistType>:
//     PolyDegree<T>
// {
//     type Error: std::error::Error + Send + Sync + 'static;
//     fn sample_hash<B: AsRef<[u8]>>(
//         &self,
//         tag: B,
//         rows: usize,
//         columns: usize,
//     ) -> Result<PolyMatrix<T, P, M>, Self::Error>;

//     fn set_key(&mut self, key: &[u8]);

//     fn expose_key(&self) -> &[u8];
// }
