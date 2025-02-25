pub mod gadget;
pub mod matrix;
pub mod sampler;
use phantom_zone_math::modulus::{Elem, ElemFrom, ElemOps};
use std::fmt::Debug;

pub type PElem<T> = Elem<T>;

pub trait PolyElemModulus {
    fn modulus_bits(&self) -> usize;
}
pub trait PolyElemOps:
    ElemOps + ElemFrom<u64> + ElemFrom<u32> + ElemFrom<u8> + ElemFrom<bool> + PolyElemModulus
{
    type Error: std::error::Error + Send + Sync + 'static;
}

// pub trait PolyBitOps: PolyElemOps + ElemTo<u64> + ElemTo<u32> + ElemTo<u8> + ElemTo<bool> {
//     fn modulus_bits(&self) -> usize {
//         1
//     }
// }

// pub trait PolyGaussOps: PolyElemOps {
//     fn gaussian_param(&self) -> f64;
// }

pub type Poly<T, P> = <P as PolyOps<T>>::Poly;

pub trait PolyDegree<T: PolyElemOps>: PolyElemModulus {
    fn degree(&self) -> usize;
}

pub trait PolyOps<T: PolyElemOps>: PolyDegree<T> {
    type Error: std::error::Error + Send + Sync + 'static;
    type Poly: Debug + Clone;
    fn coeffs(&self, poly: &Self::Poly) -> &[PElem<T>];
    fn from_coeffs(coeffs: &[PElem<T>]) -> Result<Self::Poly, Self::Error>;
    fn from_const(constant: &T) -> Result<Self::Poly, Self::Error>;
    fn zero(&self) -> Result<Self::Poly, Self::Error>;
    fn one(&self) -> Result<Self::Poly, Self::Error>;
    fn minus_one(&self) -> Result<Self::Poly, Self::Error>;
    fn add(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn neg(&self, a: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn sub(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error> {
        let minus_b = self.neg(b)?;
        self.add(a, &minus_b)
    }
    fn mul(&self, a: &Self::Poly, b: &Self::Poly) -> Result<Self::Poly, Self::Error>;
    fn extract_highest_bits(&self, poly: &Self::Poly) -> Result<Vec<bool>, Self::Error>;
}

use openfhe::{cxx::CxxVector, ffi};

fn main() {
    use std::time::Instant;
    println!("\n======EXAMPLE FOR EVALPOLY========\n");

    let mut _cc_params_ckksrns = ffi::GenParamsCKKSRNS();
    _cc_params_ckksrns.pin_mut().SetMultiplicativeDepth(6);
    _cc_params_ckksrns.pin_mut().SetScalingModSize(50);

    let _cc = ffi::DCRTPolyGenCryptoContextByParamsCKKSRNS(&_cc_params_ckksrns);
    _cc.EnableByFeature(ffi::PKESchemeFeature::PKE);
    _cc.EnableByFeature(ffi::PKESchemeFeature::KEYSWITCH);
    _cc.EnableByFeature(ffi::PKESchemeFeature::LEVELEDSHE);
    _cc.EnableByFeature(ffi::PKESchemeFeature::ADVANCEDSHE);

    let mut _input = CxxVector::<ffi::ComplexPair>::new();
    _input.pin_mut().push(ffi::ComplexPair { re: 0.5, im: 0.0 });
    _input.pin_mut().push(ffi::ComplexPair { re: 0.7, im: 0.0 });
    _input.pin_mut().push(ffi::ComplexPair { re: 0.9, im: 0.0 });
    _input.pin_mut().push(ffi::ComplexPair { re: 0.95, im: 0.0 });
    _input.pin_mut().push(ffi::ComplexPair { re: 0.93, im: 0.0 });
    let _encoded_length = _input.len();

    let mut _coefficients_1 = CxxVector::<f64>::new();
    _coefficients_1.pin_mut().push(0.15);
    _coefficients_1.pin_mut().push(0.75);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(1.25);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(1.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(1.0);
    _coefficients_1.pin_mut().push(2.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(1.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(0.0);
    _coefficients_1.pin_mut().push(1.0);

    let mut _coefficients_2 = CxxVector::<f64>::new();
    _coefficients_2.pin_mut().push(1.0);
    _coefficients_2.pin_mut().push(2.0);
    _coefficients_2.pin_mut().push(3.0);
    _coefficients_2.pin_mut().push(4.0);
    _coefficients_2.pin_mut().push(5.0);
    _coefficients_2.pin_mut().push(-1.0);
    _coefficients_2.pin_mut().push(-2.0);
    _coefficients_2.pin_mut().push(-3.0);
    _coefficients_2.pin_mut().push(-4.0);
    _coefficients_2.pin_mut().push(-5.0);
    _coefficients_2.pin_mut().push(0.1);
    _coefficients_2.pin_mut().push(0.2);
    _coefficients_2.pin_mut().push(0.3);
    _coefficients_2.pin_mut().push(0.4);
    _coefficients_2.pin_mut().push(0.5);
    _coefficients_2.pin_mut().push(-0.1);
    _coefficients_2.pin_mut().push(-0.2);
    _coefficients_2.pin_mut().push(-0.3);
    _coefficients_2.pin_mut().push(-0.4);
    _coefficients_2.pin_mut().push(-0.5);
    _coefficients_2.pin_mut().push(0.1);
    _coefficients_2.pin_mut().push(0.2);
    _coefficients_2.pin_mut().push(0.3);
    _coefficients_2.pin_mut().push(0.4);
    _coefficients_2.pin_mut().push(0.5);
    _coefficients_2.pin_mut().push(-0.1);
    _coefficients_2.pin_mut().push(-0.2);
    _coefficients_2.pin_mut().push(-0.3);
    _coefficients_2.pin_mut().push(-0.4);
    _coefficients_2.pin_mut().push(-0.5);

    let _dcrt_poly_params = ffi::DCRTPolyGenNullParams();
    let _plain_text_1 =
        _cc.MakeCKKSPackedPlaintextByVectorOfComplex(&_input, 1, 0, &_dcrt_poly_params, 0);
    let _key_pair = _cc.KeyGen();
    print!("Generating evaluation key for homomorphic multiplication...");
    _cc.EvalMultKeyGen(&_key_pair.GetPrivateKey());
    println!("Completed.\n");
    let mut _cipher_text_1 = _cc.EncryptByPublicKey(&_key_pair.GetPublicKey(), &_plain_text_1);

    let mut _start = Instant::now();
    let _result = _cc.EvalPoly(&_cipher_text_1, &_coefficients_1);
    let _time_eval_poly_1 = _start.elapsed();

    _start = Instant::now();
    let _result_2 = _cc.EvalPoly(&_cipher_text_1, &_coefficients_2);
    let _time_eval_poly_2 = _start.elapsed();

    let mut _plain_text_dec = ffi::GenNullPlainText();
    _cc.DecryptByPrivateKeyAndCiphertext(
        &_key_pair.GetPrivateKey(),
        &_result,
        _plain_text_dec.pin_mut(),
    );
    _plain_text_dec.SetLength(_encoded_length);
    let mut _plain_text_dec_2 = ffi::GenNullPlainText();
    _cc.DecryptByPrivateKeyAndCiphertext(
        &_key_pair.GetPrivateKey(),
        &_result_2,
        _plain_text_dec_2.pin_mut(),
    );
    _plain_text_dec_2.SetLength(_encoded_length);

    println!("\n Original Plaintext #1:");
    println!("{}", _plain_text_1.GetString());
    println!(
        "\n Result of evaluating a polynomial with coefficients [{} ]",
        _coefficients_1.iter().fold(String::new(), |acc, &arg| acc + " " + &arg.to_string())
    );
    println!("{}", _plain_text_dec.GetString());
    println!("\n Expected result: (0.70519107, 1.38285078, 3.97211180, 5.60215665, 4.86357575)");
    println!("\n Evaluation time: {:.0?}", _time_eval_poly_1);
    println!(
        "\n Result of evaluating a polynomial with coefficients [{} ]",
        _coefficients_2.iter().fold(String::new(), |acc, &arg| acc + " " + &arg.to_string())
    );
    println!("{}\n", _plain_text_dec_2.GetString());
    println!(
        " Expected result: (3.4515092326, 5.3752765397, 4.8993108833, 3.2495023573, 4.0485229982)"
    );
    print!("\n Evaluation time: {:.0?}\n", _time_eval_poly_2);
}
