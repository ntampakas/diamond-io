// #[allow(clippy::type_complexity)]
// pub trait PolyTrapdoorSampler<T: PolyElemOps, P: PolyOps<T>, M: PolyMatrixOps<T, P>, D:
// DistType>:     PolyDegree<T>
// {
//     type Error: std::error::Error + Send + Sync + 'static;
//     type Trapdoor;
//     fn trapdoor<R: RngCore>(
//         &self,
//         rng: &mut R,
//     ) -> Result<(PolyMatrix<T, P, M>, Self::Trapdoor), Self::Error>;

//     fn preimage<R: RngCore>(
//         &self,
//         rng: &mut R,
//         trapdoor: &Self::Trapdoor,
//         target: &PolyMatrix<T, P, M>,
//     ) -> Result<PolyMatrix<T, P, M>, Self::Error>;
// }
