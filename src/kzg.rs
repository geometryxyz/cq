use std::{iter, marker::PhantomData};

use ark_ec::{msm::VariableBaseMSM, PairingEngine};
use ark_ff::{One, PrimeField};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};

use crate::utils::x_pow_d;

pub enum DegreeBound {
    LessThan,
    Exact,
}

pub struct Kzg<E: PairingEngine> {
    _e: PhantomData<E>,
}

impl<E: PairingEngine> Kzg<E> {
    pub fn commit_g1(srs: &[E::G1Affine], poly: &DensePolynomial<E::Fr>) -> E::G1Projective {
        if srs.len() - 1 < poly.degree() {
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                srs.len()
            );
        }
        let coeff_scalars: Vec<_> = poly.coeffs.iter().map(|c| c.into_repr()).collect();
        VariableBaseMSM::multi_scalar_mul(&srs, &coeff_scalars)
    }

    pub fn commit_g2(srs: &[E::G2Affine], poly: &DensePolynomial<E::Fr>) -> E::G2Projective {
        if srs.len() - 1 < poly.degree() {
            panic!(
                "SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                poly.degree(),
                srs.len()
            );
        }
        let coeff_scalars: Vec<_> = poly.coeffs.iter().map(|c| c.into_repr()).collect();
        VariableBaseMSM::multi_scalar_mul(&srs, &coeff_scalars)
    }

    pub fn open_g1(
        srs: &[E::G1Affine],
        poly: &DensePolynomial<E::Fr>,
        challenge: E::Fr,
    ) -> (E::Fr, E::G1Affine) {
        let q = poly / &DensePolynomial::from_coefficients_slice(&[-challenge, E::Fr::one()]);
        if srs.len() - 1 < q.degree() {
            panic!(
                "Open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                srs.len()
            );
        }
        let proof = Self::commit_g1(srs, &q);
        (poly.evaluate(&challenge), proof.into())
    }

    pub fn batch_open_g1(
        srs: &[E::G1Affine],
        polys: &[DensePolynomial<E::Fr>],
        opening_challenge: E::Fr,
        separation_challenge: E::Fr,
    ) -> E::G1Affine {
        let powers_of_gamma = iter::successors(Some(separation_challenge), |p| {
            Some(p.clone() * separation_challenge)
        });

        let mut batched = polys[0].clone();
        for (p_i, gamma_pow_i) in polys.iter().skip(1).zip(powers_of_gamma) {
            batched += (gamma_pow_i, p_i);
        }

        let q = &batched
            / &DensePolynomial::from_coefficients_slice(&[-opening_challenge, E::Fr::one()]);

        if srs.len() - 1 < q.degree() {
            panic!(
                "Batch open g1: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                srs.len()
            );
        }

        Self::commit_g1(srs, &q).into()
    }

    /// Specific function for step 15
    pub fn batch_open_g1_at_zero_with_bounds(
        srs: &[E::G1Affine],
        polys: &[DensePolynomial<E::Fr>],
        poly_bounds: &[DegreeBound],
        separation_challenge: E::Fr,
        bound: usize,
    ) -> E::G1Affine {
        let n = polys.len();
        let powers_of_gamma: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| {
            Some(p.clone() * separation_challenge)
        })
        .take(2 * n)
        .collect();

        let mut batched = DensePolynomial::default();
        for (i, p_i) in polys.iter().enumerate() {
            batched += (powers_of_gamma[i], p_i);
        }

        batched = DensePolynomial::from_coefficients_slice(&batched.coeffs[1..]);

        let d = srs.len() - 1;
        let s = d - bound + 1;

        let x_exact = -x_pow_d::<E::Fr>(bound);
        let x_upper = x_pow_d::<E::Fr>(s);

        for (i, p_i) in polys.iter().enumerate() {
            let poly = match poly_bounds[i] {
                DegreeBound::LessThan => p_i.clone(),
                DegreeBound::Exact => p_i + &x_exact,
            };
            batched += (powers_of_gamma[n + i], &(&x_upper * &poly));
        }

        if srs.len() - 1 < batched.degree() {
            panic!(
                "Batch open g1 at zero with bounds: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                batched.degree(),
                srs.len()
            );
        }

        Self::commit_g1(srs, &batched).into()
    }

    pub fn open_g2(
        srs: &[E::G2Affine],
        poly: &DensePolynomial<E::Fr>,
        challenge: E::Fr,
    ) -> (E::Fr, E::G2Affine) {
        let q = poly / &DensePolynomial::from_coefficients_slice(&[-challenge, E::Fr::one()]);
        if srs.len() - 1 < q.degree() {
            panic!(
                "Open g2: SRS size to small! Can't commit to polynomial of degree {} with srs of size {}",
                q.degree(),
                srs.len()
            );
        }
        let proof = Self::commit_g2(srs, &q);
        (poly.evaluate(&challenge), proof.into())
    }
}
