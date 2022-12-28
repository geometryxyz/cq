use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::Field;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
};

use crate::kzg::Kzg;

pub fn compute_qs<E: PairingEngine>(
    t: &DensePolynomial<E::Fr>,
    domain: &GeneralEvaluationDomain<E::Fr>,
    srs_g1: &[E::G1Affine],
) -> Vec<E::G1Affine> {
    // TODO: implement this in NlogN with the algorithm of Feist and Khovratovich with matrix product

    let n = domain.size();
    assert!(t.degree() < n);
    assert_eq!(srs_g1.len(), n); // x^0, ..., x^(n-1)

    let roots: Vec<_> = domain.elements().collect();
    let ks: Vec<E::G1Affine> = roots
        .iter()
        .map(|&g_i| Kzg::<E>::open_g1(&srs_g1, t, g_i).1.into())
        .collect();

    let n_inv = domain.size_as_field_element().inverse().unwrap();
    let normalized_roots: Vec<_> = roots.iter().map(|&g_i| g_i * n_inv).collect();

    ks.iter()
        .zip(normalized_roots.iter())
        .map(|(ki, &normalizer_i)| ki.mul(normalizer_i).into())
        .collect()
}
