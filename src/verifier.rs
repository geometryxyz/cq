use std::{iter, marker::PhantomData, ops::Neg};

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::{Field, One};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

use crate::{
    data_structures::{Proof, Statement},
    error::Error,
    indexer::CommonPreprocessedInput,
    rng::FiatShamirRng,
    transcript::TranscriptOracle,
    PROTOCOL_NAME,
};

pub struct VerifierKey<E: PairingEngine> {
    pub(crate) x: E::G2Prepared,
    pub(crate) x_pow_b0_bound: E::G2Prepared,
    pub(crate) table_size: usize,
    pub(crate) witness_size: usize,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub fn new(srs_g2: &[E::G2Affine], table_size: usize, witness_size: usize) -> Self {
        Self {
            x: srs_g2[1].into(),
            x_pow_b0_bound: srs_g2[table_size - witness_size + 1].into(),
            table_size,
            witness_size,
        }
    }
}

pub struct Verifier<E: PairingEngine, FS: FiatShamirRng> {
    _e: PhantomData<E>,
    _fs: PhantomData<FS>,
}

impl<E: PairingEngine, FS: FiatShamirRng> Verifier<E, FS> {
    pub fn verify(
        vk: &VerifierKey<E>,
        common: &CommonPreprocessedInput<E>,
        statement: &Statement<E>,
        proof: &Proof<E>,
    ) -> Result<(), Error> {
        let mut transcipt = TranscriptOracle::<FS>::initialize(&PROTOCOL_NAME);

        transcipt.stream_public_input(common, statement);

        transcipt.stream_first_message(&proof.first_msg);

        let beta: E::Fr = transcipt.squeeze_challenge();

        transcipt.stream_second_message(&proof.second_msg);

        let gamma: E::Fr = transcipt.squeeze_challenge();
        let eta: E::Fr = transcipt.squeeze_challenge();

        transcipt.stream_third_message(&proof.third_msg);

        // separator for pairing batching
        let u: E::Fr = transcipt.squeeze_challenge();
        let u_powers: Vec<E::Fr> = iter::successors(Some(u), |u_pow| Some(*u_pow * u))
            .take(4)
            .collect();

        // NOTE: for easier convention, every pairing that is written on rhs of paper will be negated for usage in product of pairings

        let g_1 = E::G1Affine::prime_subgroup_generator();
        let g_2 = E::G2Affine::prime_subgroup_generator();

        let witness_domain = GeneralEvaluationDomain::<E::Fr>::new(vk.witness_size).unwrap();

        let n_table = E::Fr::from(vk.table_size as u64);
        let n = E::Fr::from(vk.witness_size as u64);

        let b0 = n_table * proof.third_msg.a_at_zero * n.inverse().unwrap();
        let b_at_gamma = proof.third_msg.b0_at_gamma * gamma + b0;
        let f_at_gamma = proof.third_msg.f_at_gamma;
        let zh_at_gamma = witness_domain.evaluate_vanishing_polynomial(gamma);

        let qb_at_gamma =
            (b_at_gamma * (f_at_gamma + beta) - E::Fr::one()) * zh_at_gamma.inverse().unwrap();

        let v = proof.third_msg.b0_at_gamma + eta * f_at_gamma + eta * eta * qb_at_gamma;
        let minus_v_g1 = g_1.mul(-v).into_affine();
        let mut c = statement.f.mul(eta) + proof.second_msg.qb_cm.mul(eta * eta);
        c.add_assign_mixed(&proof.second_msg.b0_cm);
        let c = c.into_affine();

        let l: E::G1Affine = proof
            .third_msg
            .pi_gamma
            .mul(gamma)
            .add_mixed(&(c + minus_v_g1))
            .into();
        let minus_a_at_zero = g_1.mul(-proof.third_msg.a_at_zero).into_affine();
        let a_pt = proof.second_msg.a_cm + minus_a_at_zero;

        let beta_2 = g_2.mul(beta).into_affine();

        let lhs_batched_1 = (proof.first_msg.m_cm.neg().into_projective()
            + proof.second_msg.p_cm.mul(-u_powers[0])
            + l.mul(u_powers[1])
            + a_pt.mul(u_powers[2]))
        .into_affine();
        let lhs_batched_x = (proof.third_msg.pi_gamma.mul(u_powers[1])
            + proof.third_msg.a0_cm.mul(u_powers[2]))
        .neg()
        .into_affine();

        let res = E::product_of_pairings(&[
            (lhs_batched_1.into(), g_2.into()),
            (lhs_batched_x.into(), vk.x.clone()),
            (
                proof.second_msg.b0_cm.mul(u_powers[0]).into_affine().into(),
                vk.x_pow_b0_bound.clone(),
            ),
            (proof.second_msg.qa_cm.neg().into(), common.zv_2.into()),
            (proof.second_msg.a_cm.into(), (common.t_2 + beta_2).into()),
        ]);

        if res != E::Fqk::one() {
            if cfg!(feature = "debug") {
                // check well formation of A
                {
                    let res = E::product_of_pairings(&[
                        (proof.second_msg.a_cm.into(), (common.t_2 + beta_2).into()),
                        (proof.second_msg.qa_cm.neg().into(), common.zv_2.into()),
                        (proof.first_msg.m_cm.neg().into(), g_2.into()),
                    ]);

                    if res != E::Fqk::one() {
                        return Err(Error::Pairing1Failed);
                    }
                }

                // check b0 degree
                {
                    let res = E::product_of_pairings(&[
                        (proof.second_msg.b0_cm.into(), vk.x_pow_b0_bound.clone()),
                        (proof.second_msg.p_cm.neg().into(), g_2.into()),
                    ]);

                    if res != E::Fqk::one() {
                        return Err(Error::Pairing2Failed);
                    }
                }

                // check openings at gamma
                {
                    let res = E::product_of_pairings(&[
                        (l.into(), g_2.into()),
                        (proof.third_msg.pi_gamma.neg().into(), vk.x.clone()),
                    ]);

                    if res != E::Fqk::one() {
                        return Err(Error::Pairing3Failed);
                    }
                }

                // check a opening at zero
                {
                    let res = E::product_of_pairings(&[
                        (a_pt.into(), g_2.into()),
                        (proof.third_msg.a0_cm.neg().into(), vk.x.clone()),
                    ]);

                    if res != E::Fqk::one() {
                        return Err(Error::Pairing4Failed);
                    }
                }
            } else {
                return Err(Error::BatchedPairingFailed);
            }
        }

        Ok(())
    }
}
