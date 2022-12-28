use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::{Field, One, PrimeField, FftField};
use ark_poly::{univariate::DensePolynomial, UVPolynomial};
use ark_std::rand::RngCore;
use ark_std::UniformRand;
use std::{cmp::max, iter};

/// Create srs from rng
pub fn unsafe_setup_from_rng<E: PairingEngine, R: RngCore>(
    max_power_g1: usize,
    max_power_g2: usize,
    rng: &mut R,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let tau = E::Fr::rand(rng);
    let size = max(max_power_g1 + 1, max_power_g2 + 1);
    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(p.clone() * tau))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(max_power_g1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(max_power_g2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}

/// Create srs from specific tau
pub fn unsafe_setup_from_tau<E: PairingEngine, R: RngCore>(
    max_power_g1: usize,
    max_power_g2: usize,
    tau: E::Fr,
) -> (Vec<E::G1Affine>, Vec<E::G2Affine>) {
    let size = max(max_power_g1 + 1, max_power_g2 + 1);
    let powers_of_tau: Vec<E::Fr> = iter::successors(Some(E::Fr::one()), |p| Some(p.clone() * tau))
        .take(size)
        .collect();

    let g1_gen = E::G1Affine::prime_subgroup_generator();
    let g2_gen = E::G2Affine::prime_subgroup_generator();

    let srs_g1: Vec<E::G1Affine> = powers_of_tau
        .iter()
        .take(max_power_g1 + 1)
        .map(|tp| g1_gen.mul(tp.into_repr()).into())
        .collect();

    let srs_g2: Vec<E::G2Affine> = powers_of_tau
        .iter()
        .take(max_power_g2 + 1)
        .map(|tp| g2_gen.mul(tp.into_repr()).into())
        .collect();
    (srs_g1, srs_g2)
}

// given x coords construct Li polynomials
pub fn construct_lagrange_basis<F: FftField>(evaluation_domain: &[F]) -> Vec<DensePolynomial<F>> {
    let mut bases = Vec::with_capacity(evaluation_domain.len());
    for i in 0..evaluation_domain.len() {
        let mut l_i = DensePolynomial::from_coefficients_slice(&[F::one()]);
        let x_i = evaluation_domain[i];
        for j in 0..evaluation_domain.len() {
            if j != i {
                let xi_minus_xj_inv = (x_i - evaluation_domain[j]).inverse().unwrap();
                l_i = &l_i
                    * &DensePolynomial::from_coefficients_slice(&[
                        -evaluation_domain[j] * xi_minus_xj_inv,
                        xi_minus_xj_inv,
                    ]);
            }
        }

        bases.push(l_i);
    }

    bases
}

pub fn x_pow_d<F: Field>(d: usize) -> DensePolynomial<F> {
    let mut coeffs = vec![F::zero(); d];
    coeffs.push(F::one());
    DensePolynomial::from_coefficients_slice(&coeffs)
}

pub fn is_pow_2(x: usize) -> bool {
    (x & (x - 1)) == 0
}
