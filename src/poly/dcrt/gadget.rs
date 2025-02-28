use num_bigint::BigInt;

use crate::{
    poly::{
        dcrt::{DCRTPoly, FinRing},
        gadget::PolyGadget,
        Poly, PolyMatrix, PolyParams,
    },
    utils::ceil_log2,
};

use super::{DCRTPolyMatrix, DCRTPolyParams};

pub struct DCRTPolyGadget {}

impl PolyGadget for DCRTPolyGadget {
    type Error = std::io::Error;
    type M = DCRTPolyMatrix;
    type Params = DCRTPolyParams;

    /// Gadget vector g = (2^0, 2^1, ..., 2^{log(q)-1})
    /// where g ∈ Z_q^{log(q)}
    fn gadget_vector(params: &Self::Params) -> Self::M {
        let q = params.modulus();
        let size = ceil_log2(&q) + 1;
        let mut poly_vec = Vec::with_capacity(size);
        for i in 0..size {
            let value = BigInt::from(2).pow(i.try_into().unwrap());
            let fe: FinRing = FinRing::new(value, q.clone().into());
            poly_vec.push(DCRTPoly::from_const(&params, &fe));
        }
        DCRTPolyMatrix::from_poly_vec(&params, vec![poly_vec])
    }

    fn gadget_matrix(params: &Self::Params, size: usize) -> Self::M {
        let identity = DCRTPolyMatrix::identity(&params, size, None);
        let gadget_vector = Self::gadget_vector(&params);
        identity.tensor(&gadget_vector.transpose())
    }

    fn decompose(&self) -> Result<Self::M, Self::Error> {
        todo!()
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::poly::dcrt::DCRTPolyParams;

    #[test]
    fn test_gadget_vector() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let gadget_vector = DCRTPolyGadget::gadget_vector(&params);
        assert_eq!(gadget_vector.row_size(), 1);
        assert_eq!(gadget_vector.col_size(), ceil_log2(&params.modulus()) + 1);
    }

    #[test]
    fn test_gadget_matrix() {
        let params = DCRTPolyParams::new(16, 4, 51);
        let size = 3;
        let gadget_matrix = DCRTPolyGadget::gadget_matrix(&params, size);
        assert_eq!(gadget_matrix.row_size(), size * (ceil_log2(&params.modulus()) + 1));
        assert_eq!(gadget_matrix.col_size(), size);
    }
}
