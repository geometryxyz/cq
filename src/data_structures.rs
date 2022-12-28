use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub struct Witness<F: Field> {
    pub(crate) size: usize,
    pub(crate) f: DensePolynomial<F>, 
    pub(crate) f_evals: Vec<F>
}