use std::{collections::BTreeMap, marker::PhantomData};

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};

use crate::{
    data_structures::{ProvingKey, Witness},
    error::Error,
    indexer::Index,
    kzg::Kzg,
    table::{self, Table},
    utils::{x_pow_d, construct_lagrange_basis},
};

pub struct Prover<E: PairingEngine> {
    _e: PhantomData<E>,
}

pub struct State<'a, E: PairingEngine> {
    pk: &'a ProvingKey<E>,
    index: &'a Index<E>,
    table: &'a Table<E::Fr>,
    witness: &'a Witness<E::Fr>,

    // captured in round_1
    m_evals: Option<BTreeMap<usize, E::Fr>>,
}

impl<'a, E: PairingEngine> State<'a, E> {
    pub fn new(
        pk: &'a ProvingKey<E>,
        index: &'a Index<E>,
        table: &'a Table<E::Fr>,
        witness: &'a Witness<E::Fr>,
    ) -> Self {
        Self {
            pk,
            index,
            table,
            witness,
            m_evals: None,
        }
    }
}

impl<E: PairingEngine> Prover<E> {
    pub fn prove() {}

    pub fn round_1<'a>(state: &'a mut State<E>) -> Result<E::G1Affine, Error> {
        let mut index_multiplicity_mapping = BTreeMap::<usize, E::Fr>::default();

        for fi in &state.witness.f_evals {
            let index = state.table.value_index_mapping.get(fi);
            let index = index.ok_or(Error::ValueNotInTable(format!("{}", fi)))?;
            let multiplicity = index_multiplicity_mapping
                .entry(*index)
                .or_insert(E::Fr::zero());
            *multiplicity = *multiplicity + E::Fr::one();
        }

        let mut m_cm = E::G1Affine::zero();
        for (&index, &multiplicity) in index_multiplicity_mapping.iter() {
            m_cm = state.index.ls[index]
                .mul(multiplicity)
                .add_mixed(&m_cm)
                .into();
        }

        state.m_evals = Some(index_multiplicity_mapping);
        Ok(m_cm)
    }

    pub fn round_2<'a>(
        state: &'a mut State<E>,
        beta: E::Fr,
    ) -> Result<
        (
            E::G1Affine,
            E::G1Affine,
            E::G1Affine,
            E::G1Affine,
            E::G1Affine,
        ),
        Error,
    > {
        let wtns_domain = GeneralEvaluationDomain::<E::Fr>::new(state.witness.size).unwrap();
        let m_evals = state.m_evals.as_ref().expect("m is missing from state");

        let mut a_mapping = BTreeMap::<usize, E::Fr>::default();
        let mut a_cm = E::G1Affine::zero();

        // step 2&3: computes A sparse representation and a commitment in single pass
        for (&index, &multiplicity) in m_evals.iter() {
            let a_i = multiplicity * (state.table.values[index] + beta).inverse().unwrap();
            let _ = a_mapping.insert(index, a_i); // keys are unique so overriding will never occur

            a_cm = state.index.ls[index].mul(a_i).add_mixed(&a_cm).into();
        }

        // TODO: introduce cfg=test, and ask ifcfg = test
        // sanity 
        // {
        //     let table_domain = GeneralEvaluationDomain::<E::Fr>::new(state.table.size).unwrap();
        //     let zv: DensePolynomial<_> = table_domain.vanishing_polynomial().into();

        //     let roots: Vec<_> = table_domain.elements().collect();
        //     let lagrange_basis = construct_lagrange_basis(&roots);

        //     let mut m_poly = DensePolynomial::zero();
        //     for (&index, &multiplicity) in m_evals.iter() {
        //         m_poly += (multiplicity, &lagrange_basis[index]);
        //     }

        //     let mut a_poly = DensePolynomial::zero();
        //     for (&index, &a_i) in a_mapping.iter() {
        //         a_poly += (a_i, &lagrange_basis[index]);
        //     }

        //     let mut table_poly = DensePolynomial::from_coefficients_slice(&table_domain.ifft(&state.table.values));

        //     table_poly[0] += beta;
        //     let mut num = &a_poly * &table_poly;
        //     num += (-E::Fr::one(), &m_poly);

        //     let qa = &num / &zv; 
        //     assert_eq!(num, &qa * &zv);
        // }

        // step 4: compute [QA(x)]_1
        let mut qa_cm = E::G1Affine::zero();
        for (&index, &a_i) in a_mapping.iter() {
            qa_cm = state.index.qs[index].mul(a_i).add_mixed(&qa_cm).into();
        }


        // step 5: compute B(X)
        let b_evals: Vec<_> = state
            .witness
            .f_evals
            .iter()
            .map(|&fi| (fi + beta).inverse().unwrap())
            .collect();
        let b_poly = DensePolynomial::from_coefficients_slice(&wtns_domain.ifft(&b_evals));

        // step 6: compute B0(X)
        let b0_poly = DensePolynomial::from_coefficients_slice(&b_poly.coeffs[1..]);

        // step 7: commit to B0(X)
        let b0_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &b0_poly).into();

        // step 8: compute QB(X)
        let qb_evals: Vec<_> = b_evals
            .iter()
            .zip(state.witness.f_evals.iter())
            .map(|(&bi, &fi)| bi * (fi + beta) - E::Fr::one())
            .collect();
        let qb_poly = DensePolynomial::from_coefficients_slice(&wtns_domain.ifft(&qb_evals));

        // step 9: commit to QB(X)
        let qb_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &qb_poly).into();

        // step 10: compute degree correctness check for B0
        let p_poly = &b0_poly * &x_pow_d(state.table.size - (state.witness.size + 1));
        let p_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &p_poly).into();

        Ok((a_cm, qa_cm, b0_cm, qb_cm, p_cm))
    }
}

#[cfg(test)]
mod prover_rounds_tests {
    use std::ops::Neg;

    use ark_bn254::{Bn254, Fr, G2Affine, Fq12};
    use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
    use ark_ff::{UniformRand, One};
    use ark_std::{rand::rngs::StdRng, test_rng};

    use crate::{
        data_structures::{ProvingKey, Witness},
        indexer::Index,
        table::Table,
        utils::{to_field, unsafe_setup_from_rng},
    };

    use super::{Prover, State};

    // TODO: create prepare function for shared data generation

    #[test]
    fn test_round_1() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1, srs_g2 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &pk.srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let mut state = State::new(&pk, &index, &table, &witness);

        let res = Prover::round_1(&mut state);
        assert!(res.is_ok());

        let keys = vec![1, 3, 4, 7];
        let supp_m: Vec<usize> = state.m_evals.as_ref().unwrap().keys().map(|&i| i).collect();
        assert_eq!(keys, supp_m);

        let multiplicities = vec![Fr::one(), Fr::one(), Fr::one(), Fr::one()];
        let m_values: Vec<Fr> = state.m_evals.as_ref().unwrap().values().map(|&mi| mi).collect();
        assert_eq!(multiplicities, m_values);
    }

    #[test]
    pub fn test_round_2() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1, srs_g2 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &pk.srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let mut state = State::new(&pk, &index, &table, &witness);

        let m_cm = Prover::round_1(&mut state).unwrap();

        let beta = Fr::rand(&mut rng);
        let (a_cm, qa_cm, b0_cm, qb_cm, p_cm) = Prover::round_2(&mut state, beta).unwrap();

        // check well formation of A
        {
            let g_2 = G2Affine::prime_subgroup_generator();
            let beta_2 = g_2.mul(beta).into_affine();
            let rhs_0 = index.t_2 + beta_2;

            let res = Bn254::product_of_pairings(&[
                (a_cm.neg().into(), rhs_0.into()), 
                (qa_cm.into(), index.zv_2.into()), 
                (m_cm.into(), g_2.into())
            ]);

            assert_eq!(res, Fq12::one());
        }

        // check b0 degree 
        {
            let lhs_0 = b0_cm;
            let rhs_0 = pk.srs_g2[table.size - witness.size - 1];

            let lhs_1 = p_cm; 
            let rhs_1 = G2Affine::prime_subgroup_generator();
            assert_eq!(
                Bn254::pairing(lhs_0, rhs_0),
                Bn254::pairing(lhs_1, rhs_1),
            )
        }
    }
}
