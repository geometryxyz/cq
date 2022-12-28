use std::{marker::PhantomData, iter, ops::Neg};

use ark_ec::{PairingEngine, AffineCurve, ProjectiveCurve};
use ark_ff::{Field, One};
use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};

use crate::{rng::FiatShamirRng, data_structures::{Statement, Proof}, transcript::TranscriptOracle, PROTOCOL_NAME, indexer::CommonPreprocessedInput, error::Error};

pub struct VerifierKey<E: PairingEngine> {
    pub(crate) x: E::G2Affine,
    pub(crate) x_pow_b0_bound: E::G2Affine,
    pub(crate) table_size: usize, 
    pub(crate) witness_size: usize
}

impl<E: PairingEngine> VerifierKey<E> {
    pub fn new(srs_g2: &[E::G2Affine], table_size: usize, witness_size: usize) -> Self {
        Self {
            x: srs_g2[1], 
            x_pow_b0_bound: srs_g2[table_size - witness_size - 1], 
            table_size, 
            witness_size
        }
    }
}

pub struct Verifier<E: PairingEngine, FS: FiatShamirRng> {
    _e: PhantomData<E>, 
    _fs: PhantomData<FS>
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
        let u_powers: Vec<E::Fr> = iter::successors(Some(u), |u_pow| Some(*u_pow * u)).take(4).collect();

        // NOTE: for easier convention, every pairing that is written on rhs of paper will be negated for usage in product of pairings

        let g_1 = E::G1Affine::prime_subgroup_generator();
        let g_2 = E::G2Affine::prime_subgroup_generator();

        let witness_domain = GeneralEvaluationDomain::<E::Fr>::new(vk.witness_size).unwrap();

        let N = E::Fr::from(vk.table_size as u64); 
        let n = E::Fr::from(vk.witness_size as u64);

        let b0 = N * proof.third_msg.a_at_zero * n.inverse().unwrap();
        let b_at_gamma = proof.third_msg.b0_at_gamma * gamma + b0;
        let f_at_gamma = proof.third_msg.f_at_gamma; 
        let zh_at_gamma = witness_domain.evaluate_vanishing_polynomial(gamma);

        let qb_at_gamma = (b_at_gamma * (f_at_gamma + beta) - E::Fr::one()) * zh_at_gamma.inverse().unwrap();

        let v = proof.third_msg.b0_at_gamma + eta * f_at_gamma + eta * eta * qb_at_gamma;
        let minus_v_g1 = g_1.mul(-v).into_affine();
        let mut c = statement.f.mul(eta) + proof.second_msg.qb_cm.mul(eta * eta);
        c.add_assign_mixed(&proof.second_msg.b0_cm);
        let c = c.into_affine();

        let l: E::G1Affine = proof.third_msg.pi_gamma.mul(gamma).add_mixed(&(c + minus_v_g1)).into();
        let minus_a_at_zero = g_1.mul(-proof.third_msg.a_at_zero).into_affine();
        let a_pt = proof.third_msg.a0_cm + minus_a_at_zero;

        let beta_2 = g_2.mul(beta).into_affine();

        // pairing 1
        {
            let g_2 = E::G2Affine::prime_subgroup_generator();
            let beta_2 = g_2.mul(beta).into_affine();
            let rhs_0 = common.t_2 + beta_2;

            let res = E::product_of_pairings(&[
                (proof.second_msg.a_cm.neg().into(), rhs_0.into()), 
                (proof.second_msg.qa_cm.into(), common.zv_2.into()), 
                (proof.first_msg.m_cm.into(), g_2.into())
            ]);

            assert_eq!(res, E::Fqk::one());
        }

        // pairing 2 
        {
            let lhs_0 = proof.second_msg.b0_cm;

            let lhs_1 = proof.second_msg.p_cm; 
            let rhs_1 = E::G2Affine::prime_subgroup_generator();
            assert_eq!(
                E::pairing(lhs_0, vk.x_pow_b0_bound),
                E::pairing(lhs_1, rhs_1),
            )
        }

        // pairing 3 
        {
            // let N = Fr::from(table.size as u64); 
            // let n_inv = Fr::from(witness.size as u64).inverse().unwrap();

            // let b0 = N * a_at_zero * n_inv;
            // let b_at_gamma = b0 + gamma * b0_at_gamma;

            // let wtns_domain = GeneralEvaluationDomain::<Fr>::new(witness.size).unwrap();
            // let zh_at_gamma_inv = wtns_domain.evaluate_vanishing_polynomial(gamma).inverse().unwrap();

            // let qb_at_gamma = (b_at_gamma * (f_at_gamma + beta) - Fr::one()) * zh_at_gamma_inv;

            // let v = b0_at_gamma + eta * f_at_gamma + eta * eta * qb_at_gamma;
            // let mut c = statement.f.mul(eta) + qb_cm.mul(eta * eta);
            // c.add_assign_mixed(&b0_cm);
            // let c = c.into_affine();

            let g_2 = E::G2Affine::prime_subgroup_generator();
            let minus_v_g1 = E::G1Affine::prime_subgroup_generator().mul(-v).into_affine();

            let lhs: E::G1Affine = proof.third_msg.pi_gamma.mul(gamma).add_mixed(&(c + minus_v_g1)).into();
            let p1 = E::pairing(lhs, g_2);
            let p2 = E::pairing(proof.third_msg.pi_gamma, vk.x);
            assert_eq!(p1, p2);

        }

        // pairing 4
        {
            let g_2 = E::G2Affine::prime_subgroup_generator();

            let lhs = E::G1Affine::prime_subgroup_generator().mul(proof.third_msg.a_at_zero).neg().add_mixed(&proof.second_msg.a_cm);

            let p1 = E::pairing(lhs, g_2);
            let p2 = E::pairing(proof.third_msg.a0_cm, vk.x);
            assert_eq!(p1, p2);
        }


        // batched lhs of pairing that is of form e(*, [1]_2)
        // let mut lhs_1_batched = proof.second_msg.p_cm.mul(-u_powers[0]) + l.mul(u_powers[1]) + a_pt.mul(u_powers[2]);
        // lhs_1_batched.add_assign_mixed(&proof.first_msg.m_cm.neg());
        // let lhs_1_batched = lhs_1_batched.into_affine();

        // //tmp 
        // let lhs_1_batched = l.mul(u_powers[1]) + a_pt.mul(u_powers[2]);
        // let lhs_1_batched = lhs_1_batched.into_affine();

        // // batched lhs of pairing that is of form e(*, [x]_2) such that separation powers match above defined lhs_1
        // let lhs_x_batched = proof.third_msg.pi_gamma.mul(-u_powers[1]) + proof.third_msg.a0_cm.mul(-u_powers[2]);

        // let res = E::product_of_pairings(&[
        //     (lhs_1_batched.into(), g_2.into()),
        //     (lhs_x_batched.into_affine().into(), vk.x.into()),
        //     // (proof.second_msg.a_cm.into(), (common.zv_2 + beta_2).into()), 
        //     // (proof.second_msg.qa_cm.neg().into(), common.zv_2.into()), 
        //     // (proof.second_msg.b0_cm.mul(u_powers[0]).into_affine().into(), vk.x_pow_b0_bound.into())
        // ]);

        // if res != E::Fqk::one() {
        //     return Err(Error::BatchedPairingFailed);
        // }

        Ok(())

    }
}