use ark_ec::PairingEngine;
use ark_ff::FftField;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};

use crate::error::Error;

pub struct ProvingKey<E: PairingEngine> {
    pub(crate) srs_g1: Vec<E::G1Affine>,
    pub(crate) srs_g2: Vec<E::G2Affine>,
}

pub struct Witness<F: FftField> {
    pub(crate) size: usize,
    pub(crate) f: DensePolynomial<F>,
    pub(crate) f_evals: Vec<F>,
}

impl<F: FftField> Witness<F> {
    pub fn new(values: &Vec<F>) -> Result<Self, Error> {
        if !values.len().is_power_of_two() {
            return Err(Error::WitnessSizeNotPow2(values.len()));
        }

        let domain = GeneralEvaluationDomain::<F>::new(values.len()).unwrap();
        let f = DensePolynomial::from_coefficients_slice(&domain.ifft(&values));

        Ok(Self {
            size: values.len(),
            f,
            f_evals: values.clone(),
        })
    }
}
