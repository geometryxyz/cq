use std::iter;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
};

use crate::{kzg::Kzg, utils::is_pow_2};

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
) -> Vec<E::G1Affine> {
    // TODO: implement this in NlogN with the algorithm of Feist and Khovratovich https://alinush.github.io/2021/06/17/Feist-Khovratovich-technique-for-computing-KZG-proofs-fast.html#mjx-eqn-eq%3Api-dft

    let n = domain.size();
    assert!(t.degree() < n);
    assert_eq!(srs_g1.len(), n); // x^0, ..., x^(n-1)

    let roots: Vec<_> = domain.elements().collect();
    let ks: Vec<E::G1Affine> = roots
        .iter()
        .map(|&g_i| Kzg::<E>::open_g1(srs_g1, t, g_i).1)
        .collect();

    let n_inv = domain.size_as_field_element().inverse().unwrap();
    let normalized_roots: Vec<_> = roots.iter().map(|&g_i| g_i * n_inv).collect();

    ks.iter()
        .zip(normalized_roots.iter())
        .map(|(ki, &normalizer_i)| ki.mul(normalizer_i).into())
        .collect()
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
