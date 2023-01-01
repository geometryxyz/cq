use std::iter;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
};
use fk::UpperToeplitz;

use crate::utils::is_pow_2;

pub fn compute_lagrange_basis_commitments<C: AffineCurve>(tau_powers: &[C]) -> Vec<C> {
    let n = tau_powers.len();
    assert!(is_pow_2(n));

    let domain = GeneralEvaluationDomain::<C::ScalarField>::new(n).unwrap();
    let n_inv = domain
        .size_as_field_element()
        .inverse()
        .unwrap()
        .into_repr();

    let tau_projective: Vec<C::Projective> = tau_powers
        .iter()
        .map(|tau_pow_i| tau_pow_i.into_projective())
        .collect();
    let p_evals: Vec<C::Projective> = domain.fft(&tau_projective);
    let p_evals_reversed: Vec<C::Projective> = iter::once(p_evals[0])
        .chain(p_evals.into_iter().skip(1).rev())
        .collect();

    let mut ls: Vec<C::Projective> = p_evals_reversed
        .into_iter()
        .map(|pi| pi.mul(n_inv))
        .collect();
    C::Projective::batch_normalization(&mut ls);
    ls.iter().map(|li| li.into_affine()).collect()
}

pub fn compute_qs<E: PairingEngine>(
    t: &DensePolynomial<E::Fr>,
    domain: &GeneralEvaluationDomain<E::Fr>,
    srs_g1: &[E::G1Affine],
)  -> Vec<E::G1Affine> {
    /* 
        - N (table size) is always pow2
        - Toeplitz multiplication will happen in 2 * N, so appending zero commitments on hs is not needed
    */

    let toeplitz = UpperToeplitz::from_poly(t);

    let mut srs_proj: Vec<E::G1Projective> = srs_g1
    .iter()
    .map(|t| t.into_projective())
    .collect();
    srs_proj.reverse();

    let h_commitments: Vec<E::G1Projective> = toeplitz.mul_by_vec(&srs_proj);
    assert_eq!(h_commitments.len(), 2 * domain.size());

    let ks: Vec<_ > = domain.fft(&toeplitz.mul_by_vec(&srs_proj)[..domain.size()]);

    let n_inv = domain.size_as_field_element().inverse().unwrap();
    let normalized_roots = domain.elements().map(|g_i| g_i * n_inv);

    let mut qs: Vec<E::G1Projective> = ks.iter()
        .zip(normalized_roots)
        .map(|(ki, normalizer_i)| ki.mul(normalizer_i.into_repr()))
        .collect();

    E::G1Projective::batch_normalization(&mut qs);
    qs.iter().map(|qi| qi.into_affine()).collect()
}

#[cfg(test)]
mod test_tools {
    use ark_bn254::{Bn254, Fr, G1Affine};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::{rand::rngs::StdRng, test_rng};

    use crate::{
        kzg::Kzg,
        utils::{construct_lagrange_basis, unsafe_setup_from_rng},
    };

    use super::compute_lagrange_basis_commitments;

    #[test]
    fn test_li_commitments() {
        let n = 2usize;
        let mut rng = test_rng();

        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let (srs_g1, _) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, 0, &mut rng);

        let roots: Vec<_> = domain.elements().collect();
        let lagrange_basis = construct_lagrange_basis(&roots);
        let lagrange_basis_1_slow: Vec<G1Affine> = lagrange_basis
            .iter()
            .map(|li| Kzg::<Bn254>::commit_g1(&srs_g1, li).into())
            .collect();

        let lagrange_basis_1_fast = compute_lagrange_basis_commitments(&srs_g1);
        assert_eq!(lagrange_basis_1_slow, lagrange_basis_1_fast);
    }
}
