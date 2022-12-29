use std::{collections::BTreeMap, marker::PhantomData};

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One, ToBytes, Zero};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};

use crate::{
    data_structures::{Proof, ProvingKey, Statement, Witness},
    error::Error,
    indexer::Index,
    kzg::Kzg,
    rng::FiatShamirRng,
    table::Table,
    transcript::TranscriptOracle,
    utils::x_pow_d,
    PROTOCOL_NAME,
};

pub struct Prover<E: PairingEngine, FS: FiatShamirRng> {
    _e: PhantomData<E>,
    _fs: PhantomData<FS>,
}

pub struct State<'a, E: PairingEngine> {
    pk: &'a ProvingKey<E>,
    index: &'a Index<E>,
    table: &'a Table<E::Fr>,
    witness: &'a Witness<E::Fr>,

    // captured in round_1
    m_sparse: Option<BTreeMap<usize, E::Fr>>,

    // captured in round_2
    b0: Option<DensePolynomial<E::Fr>>,
    qb: Option<DensePolynomial<E::Fr>>,
    a_sparse: Option<BTreeMap<usize, E::Fr>>,
    a_at_zero: Option<E::Fr>,
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

            m_sparse: None,

            b0: None,
            qb: None,
            a_sparse: None,
            a_at_zero: None,
        }
    }
}

pub struct ProverFirstMessage<E: PairingEngine> {
    pub(crate) m_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for ProverFirstMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.m_cm.write(&mut w)
    }
}

pub struct ProverSecondMessage<E: PairingEngine> {
    pub(crate) a_cm: E::G1Affine,
    pub(crate) qa_cm: E::G1Affine,
    pub(crate) b0_cm: E::G1Affine,
    pub(crate) qb_cm: E::G1Affine,
    pub(crate) p_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for ProverSecondMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.a_cm.write(&mut w)?;
        self.qa_cm.write(&mut w)?;
        self.b0_cm.write(&mut w)?;
        self.qb_cm.write(&mut w)?;
        self.p_cm.write(&mut w)
    }
}

pub struct ProverThirdMessage<E: PairingEngine> {
    pub(crate) b0_at_gamma: E::Fr,
    pub(crate) f_at_gamma: E::Fr,
    pub(crate) a_at_zero: E::Fr,
    pub(crate) pi_gamma: E::G1Affine,
    pub(crate) a0_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for ProverThirdMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.b0_at_gamma.write(&mut w)?;
        self.f_at_gamma.write(&mut w)?;
        self.a_at_zero.write(&mut w)?;
        self.pi_gamma.write(&mut w)?;
        self.a0_cm.write(&mut w)
    }
}

impl<E: PairingEngine, FS: FiatShamirRng> Prover<E, FS> {
    pub fn prove<'a>(
        pk: &'a ProvingKey<E>,
        index: &'a Index<E>,
        table: &'a Table<E::Fr>,
        witness: &'a Witness<E::Fr>,
        statement: &Statement<E>,
    ) -> Result<Proof<E>, Error> {
        let mut state = State::new(pk, index, table, witness);
        let mut transcipt = TranscriptOracle::<FS>::initialize(&PROTOCOL_NAME);

        transcipt.stream_public_input(&index.common, statement);

        let first_msg = Self::round_1(&mut state)?;
        transcipt.stream_first_message(&first_msg);

        let beta: E::Fr = transcipt.squeeze_challenge();

        let second_msg = Self::round_2(&mut state, beta)?;
        transcipt.stream_second_message(&second_msg);

        let gamma: E::Fr = transcipt.squeeze_challenge();
        let eta: E::Fr = transcipt.squeeze_challenge();

        let third_msg = Self::round_3(&mut state, gamma, eta)?;

        Ok(Proof {
            first_msg,
            second_msg,
            third_msg,
        })
    }

    pub fn round_1(state: &mut State<E>) -> Result<ProverFirstMessage<E>, Error> {
        let mut index_multiplicity_mapping = BTreeMap::<usize, E::Fr>::default();

        for fi in &state.witness.f_evals {
            let index = state.table.value_index_mapping.get(fi);
            let err_str = format!("{}", fi);
            let index = index.ok_or(Error::ValueNotInTable(err_str))?;
            let zero = E::Fr::zero();
            let multiplicity = index_multiplicity_mapping.entry(*index).or_insert(zero);
            *multiplicity += E::Fr::one();
        }

        let mut m_cm = E::G1Affine::zero();
        for (&index, &multiplicity) in index_multiplicity_mapping.iter() {
            m_cm = state.index.ls[index]
                .mul(multiplicity)
                .add_mixed(&m_cm)
                .into();
        }

        state.m_sparse = Some(index_multiplicity_mapping);
        Ok(ProverFirstMessage { m_cm })
    }

    pub fn round_2<'a>(
        state: &'a mut State<E>,
        beta: E::Fr,
    ) -> Result<ProverSecondMessage<E>, Error> {
        let wtns_domain = GeneralEvaluationDomain::<E::Fr>::new(state.witness.size).unwrap();
        let m_sparse = state.m_sparse.as_ref().expect("m is missing from state");

        let mut a_sparse = BTreeMap::<usize, E::Fr>::default();
        let mut a_cm = E::G1Affine::zero();
        let mut qa_cm = E::G1Affine::zero();

        // step 2&3&4: computes A sparse representation, a commitment and qa commitment in single pass
        for (&index, &multiplicity) in m_sparse.iter() {
            let a_i = multiplicity * (state.table.values[index] + beta).inverse().unwrap();
            let _ = a_sparse.insert(index, a_i); // keys are unique so overriding will never occur

            a_cm = state.index.ls[index].mul(a_i).add_mixed(&a_cm).into();
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
        let b_coset_evals = wtns_domain.coset_fft(&b_poly);
        let f_coset_evals = wtns_domain.coset_fft(&state.witness.f);
        let mut qb_evals: Vec<_> = b_coset_evals
            .iter()
            .zip(f_coset_evals.iter())
            .map(|(&bi, &fi)| bi * (fi + beta) - E::Fr::one())
            .collect();
        wtns_domain.divide_by_vanishing_poly_on_coset_in_place(&mut qb_evals);
        let qb_poly = DensePolynomial::from_coefficients_slice(&wtns_domain.coset_ifft(&qb_evals));

        // step 9: commit to QB(X)
        let qb_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &qb_poly).into();

        // step 10: compute degree correctness check for B0
        let p_poly = &b0_poly * &x_pow_d(state.table.size - (state.witness.size + 1));
        let p_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &p_poly).into();

        state.b0 = Some(b0_poly);
        state.qb = Some(qb_poly);

        let a_at_zero = {
            let b_at_zero = b_poly.evaluate(&E::Fr::zero());
            let n = E::Fr::from(state.witness.size as u64);

            let n_table_inv = E::Fr::from(state.table.size as u64).inverse().unwrap();

            n * b_at_zero * n_table_inv
        };

        state.a_at_zero = Some(a_at_zero);
        state.a_sparse = Some(a_sparse);

        Ok(ProverSecondMessage {
            a_cm,
            qa_cm,
            b0_cm,
            qb_cm,
            p_cm,
        })
    }

    pub fn round_3<'a>(
        state: &'a mut State<E>,
        gamma: E::Fr,
        eta: E::Fr,
    ) -> Result<ProverThirdMessage<E>, Error> {
        let b0 = state.b0.as_ref().expect("b0 is missing from state");
        let qb = state.qb.as_ref().expect("qb is missing from state");
        let a_sparse = state.a_sparse.as_ref().expect("a missing from state");
        let a_at_zero = state.a_at_zero.expect("a at 0 missing from the state");

        // step 2: compute openings of b0 and f
        let b0_at_gamma = b0.evaluate(&gamma);
        let f_at_gamma = state.witness.f.evaluate(&gamma);

        // step 3: compute [A0(X)]_1
        let mut a0_cm = E::G1Affine::zero();
        for (&index, &a_i) in a_sparse.iter() {
            a0_cm = state.index.ls_at_0[index].mul(a_i).add_mixed(&a0_cm).into();
        }

        // step 6: compute openings proof
        let pi_gamma: E::G1Affine = Kzg::<E>::batch_open_g1(
            &state.pk.srs_g1,
            &[b0.clone(), state.witness.f.clone(), qb.clone()],
            gamma,
            eta,
        )
        .into();

        Ok(ProverThirdMessage {
            b0_at_gamma,
            f_at_gamma,
            a_at_zero,
            pi_gamma,
            a0_cm,
        })
    }
}

#[cfg(test)]
mod prover_rounds_tests {
    use std::ops::Neg;

    use ark_bn254::{Bn254, Fq12, Fr, G1Affine, G2Affine};
    use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
    use ark_ff::{Field, One, UniformRand};
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
    use ark_std::{rand::rngs::StdRng, test_rng};
    use rand_chacha::ChaChaRng;
    use sha3::Keccak256;

    use crate::{
        data_structures::{ProvingKey, Statement, Witness},
        indexer::Index,
        kzg::Kzg,
        rng::SimpleHashFiatShamirRng,
        table::Table,
        utils::{to_field, unsafe_setup_from_rng},
    };

    use super::{Prover, ProverSecondMessage, ProverThirdMessage, State};

    type FS = SimpleHashFiatShamirRng<Keccak256, ChaChaRng>;

    // TODO: create prepare function for shared data generation

    #[test]
    fn test_full_proof() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let statement = Statement::<Bn254> {
            f: Kzg::<Bn254>::commit_g1(&pk.srs_g1, &witness.f).into(),
        };

        let _ = Prover::<Bn254, FS>::prove(&pk, &index, &table, &witness, &statement).unwrap();
    }

    #[test]
    fn test_round_1() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let mut state = State::new(&pk, &index, &table, &witness);

        let res = Prover::<Bn254, FS>::round_1(&mut state);
        assert!(res.is_ok());

        let keys = vec![1, 3, 4, 7];
        let supp_m: Vec<usize> = state
            .m_sparse
            .as_ref()
            .unwrap()
            .keys()
            .map(|&i| i)
            .collect();
        assert_eq!(keys, supp_m);

        let multiplicities = vec![Fr::one(), Fr::one(), Fr::one(), Fr::one()];
        let m_values: Vec<Fr> = state
            .m_sparse
            .as_ref()
            .unwrap()
            .values()
            .map(|&mi| mi)
            .collect();
        assert_eq!(multiplicities, m_values);
    }

    #[test]
    pub fn test_round_2() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let mut state = State::new(&pk, &index, &table, &witness);

        let m_cm = Prover::<Bn254, FS>::round_1(&mut state).unwrap().m_cm;

        let beta = Fr::rand(&mut rng);
        let second_msg = Prover::<Bn254, FS>::round_2(&mut state, beta).unwrap();
        let ProverSecondMessage {
            a_cm,
            qa_cm,
            b0_cm,
            qb_cm: _,
            p_cm,
        } = second_msg;

        // check well formation of A
        {
            let g_2 = G2Affine::prime_subgroup_generator();
            let beta_2 = g_2.mul(beta).into_affine();
            let rhs_0 = index.common.t_2 + beta_2;

            let res = Bn254::product_of_pairings(&[
                (a_cm.neg().into(), rhs_0.into()),
                (qa_cm.into(), index.common.zv_2.into()),
                (m_cm.into(), g_2.into()),
            ]);

            assert_eq!(res, Fq12::one());
        }

        // check b0 degree
        {
            let lhs_0 = b0_cm;
            let rhs_0 = srs_g2[table.size - witness.size - 1];

            let lhs_1 = p_cm;
            let rhs_1 = G2Affine::prime_subgroup_generator();
            assert_eq!(Bn254::pairing(lhs_0, rhs_0), Bn254::pairing(lhs_1, rhs_1),)
        }
    }

    #[test]
    fn test_round_3() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let statement = Statement::<Bn254> {
            f: Kzg::<Bn254>::commit_g1(&pk.srs_g1, &witness.f).into(),
        };

        let mut state = State::new(&pk, &index, &table, &witness);

        let _ = Prover::<Bn254, FS>::round_1(&mut state).unwrap();

        let beta = Fr::rand(&mut rng);
        let second_msg = Prover::<Bn254, FS>::round_2(&mut state, beta).unwrap();
        let ProverSecondMessage {
            a_cm,
            qa_cm: _,
            b0_cm,
            qb_cm,
            p_cm: _,
        } = second_msg;

        let gamma = Fr::rand(&mut rng);
        let eta = Fr::rand(&mut rng);

        let third_msg = Prover::<Bn254, FS>::round_3(&mut state, gamma, eta).unwrap();
        let ProverThirdMessage {
            b0_at_gamma,
            f_at_gamma,
            a_at_zero,
            pi_gamma,
            a0_cm,
        } = third_msg;

        // verifier part
        {
            let n_table = Fr::from(table.size as u64);
            let n_inv = Fr::from(witness.size as u64).inverse().unwrap();

            let b0 = n_table * a_at_zero * n_inv;
            let b_at_gamma = b0 + gamma * b0_at_gamma;

            let wtns_domain = GeneralEvaluationDomain::<Fr>::new(witness.size).unwrap();
            let zh_at_gamma_inv = wtns_domain
                .evaluate_vanishing_polynomial(gamma)
                .inverse()
                .unwrap();

            let qb_at_gamma = (b_at_gamma * (f_at_gamma + beta) - Fr::one()) * zh_at_gamma_inv;

            let v = b0_at_gamma + eta * f_at_gamma + eta * eta * qb_at_gamma;
            let mut c = statement.f.mul(eta) + qb_cm.mul(eta * eta);
            c.add_assign_mixed(&b0_cm);
            let c = c.into_affine();

            let g_2 = G2Affine::prime_subgroup_generator();
            let minus_v_g1 = G1Affine::prime_subgroup_generator()
                .mul(v)
                .into_affine()
                .neg();

            let lhs: G1Affine = pi_gamma.mul(gamma).add_mixed(&(c + minus_v_g1)).into();
            let p1 = Bn254::pairing(lhs, g_2);
            let p2 = Bn254::pairing(pi_gamma, srs_g2[1]);
            assert_eq!(p1, p2);
        }

        // check a correctness
        {
            let g_2 = G2Affine::prime_subgroup_generator();

            let lhs = G1Affine::prime_subgroup_generator()
                .mul(a_at_zero)
                .neg()
                .add_mixed(&a_cm);

            let p1 = Bn254::pairing(lhs, g_2);
            let p2 = Bn254::pairing(a0_cm, srs_g2[1]);
            assert_eq!(p1, p2);
        }
    }
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
