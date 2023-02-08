use std::{iter, marker::PhantomData, ops::Neg};

use ark_ec::{msm::VariableBaseMSM, AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{One, PrimeField};
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_std::Zero;

/// Minimal KZG functionalities needed for cq
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
        VariableBaseMSM::multi_scalar_mul(srs, &coeff_scalars)
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
        VariableBaseMSM::multi_scalar_mul(srs, &coeff_scalars)
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
            Some(*p * separation_challenge)
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
}

pub(crate) struct KzgEvaluationProof<E: PairingEngine> {
    /// Commitment
    p: E::G1Affine,
    /// Opening proof
    q: E::G1Affine,
    opening_challenge: E::Fr,
    opening: E::Fr,
}

impl<E: PairingEngine> KzgEvaluationProof<E> {
    pub fn new(
        commitment: E::G1Affine,
        proof: E::G1Affine,
        opening_challenge: E::Fr,
        opening: E::Fr,
    ) -> Self {
        Self {
            p: commitment,
            q: proof,
            opening_challenge,
            opening,
        }
    }

    fn prepare(&self) -> E::G1Projective {
        let g1_gen = E::G1Affine::prime_subgroup_generator();
        self.q.mul(self.opening_challenge).add_mixed(&(self.p)) + g1_gen.mul(self.opening).neg()
    }
}

/*
    Modified form of KZG:
    e(C - y, [1]) = e(proof, [s - x])
    e(C - y, [1]) = e(proof, s) * e(proof, -x)
    e(C - y, [1]) = e(proof, s) * e(-proof, x)
    e(C - y, [1]) * e(proof, x) = e(proof, s)
    e(C - y, [1]) * e(x * proof, [1]) = e(proof, s)
    e(C - y + x * proof, [1]) = e(proof, s)
    e(C - y + x * proof, [1]) * e(-proof, s) = g_fk^0 = 1

    Now in order to batch:
    e(C_1 - y_1 + x_1 * proof_1 + r(C_2 - y_2 + x_2 * proof_2), [1]) * e(-(proof_1 + r * proof_2), s) = g_fk^0 = 1
*/
pub(crate) fn batch_pairings<E: PairingEngine>(
    evaluation_proofs: &[KzgEvaluationProof<E>],
    r: E::Fr,
) -> (E::G1Affine, E::G1Affine) {
    let powers_of_r = iter::successors(Some(E::Fr::one()), |r_pow| Some(r_pow.clone() * r));

    let mut lhs = E::G1Projective::zero();
    let mut rhs = E::G1Projective::zero();

    for (eval_proof, r_pow) in evaluation_proofs.iter().zip(powers_of_r) {
        let r_pow_rep = r_pow.into_repr();

        let lhs_i = eval_proof.prepare();
        lhs = lhs + lhs_i.mul(r_pow_rep.clone());

        rhs = rhs + eval_proof.q.mul(r_pow_rep.clone());
    }

    (lhs.into_affine(), rhs.into_affine())
}
