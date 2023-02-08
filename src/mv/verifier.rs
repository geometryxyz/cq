use std::{marker::PhantomData, ops::Neg};

use ark_ec::{AffineCurve, PairingEngine};
use ark_ff::ToBytes;
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_std::One;

use crate::{
    error::Error,
    kzg::{batch_pairings, KzgEvaluationProof},
    rng::FiatShamirRng,
};

use super::{
    data_structures::Statement, prover::Proof, transcript::TranscriptOracle, MV_PROTOCOL_NAME,
};

#[derive(Clone, Copy)]
pub struct VerifierKey<E: PairingEngine> {
    pub(crate) x: E::G2Prepared,
    pub(crate) table_cm: E::G1Affine,
    pub(crate) table_size: usize,
    pub(crate) g2_gen: E::G2Prepared,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub fn new(srs_g2: &[E::G2Affine], table_size: usize, table_cm: E::G1Affine) -> Self {
        Self {
            x: srs_g2[1].into(),
            table_size,
            table_cm,
            g2_gen: E::G2Affine::prime_subgroup_generator().into(),
        }
    }
}

impl<E: PairingEngine> ToBytes for VerifierKey<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.x.write(&mut w)?;
        E::Fr::from(self.table_size as u64).write(&mut w)?;
        self.table_cm.write(&mut w)?;
        self.g2_gen.write(&mut w)
    }
}

pub struct Verifier<E: PairingEngine, FS: FiatShamirRng>(PhantomData<(E, FS)>);

impl<E: PairingEngine, FS: FiatShamirRng> Verifier<E, FS> {
    pub fn verify(
        vk: &VerifierKey<E>,
        statement: &Statement<E>,
        proof: &Proof<E>,
    ) -> Result<(), Error> {
        let mut transcript = TranscriptOracle::<FS>::initialize(&MV_PROTOCOL_NAME);

        transcript.stream_public_input(vk, statement);
        transcript.stream_first_message(&proof.first_msg);

        let beta: E::Fr = transcript.squeeze_challenge();

        transcript.stream_second_message(&proof.second_msg);

        let x: E::Fr = transcript.squeeze_challenge();

        transcript.stream_third_message(&proof.third_msg);
        let r: E::Fr = transcript.squeeze_challenge();

        let domain = GeneralEvaluationDomain::<E::Fr>::new(vk.table_size).unwrap();

        // Check grand sum identity
        // phi(Xg) - phi(X) = 1/(beta + f) - m/(beta + t)
        // (beta + t) * (beta + f) * (phi(Xg) - phi(X)) - ((beta + t) - m*(beta + f)) = Q(X) * ZH
        {
            let zh_at_x = domain.evaluate_vanishing_polynomial(x);

            let (t_eval, f_eval, m_eval, phi_eval, phi_shifted_eval, q_eval) = (
                proof.third_msg.t_opening.0,
                proof.third_msg.f_opening.0,
                proof.third_msg.m_opening.0,
                proof.third_msg.phi_opening.0,
                proof.third_msg.phi_shifted_opening.0,
                proof.third_msg.q_opening.0,
            );

            let lhs = (beta + t_eval) * (beta + f_eval) * (phi_shifted_eval - phi_eval);
            let rhs = (beta + t_eval) - (m_eval * (beta + f_eval));

            if lhs - rhs != q_eval * zh_at_x {
                return Err(Error::MVGrandSumFailed);
            }
        }

        // Check openings
        let f_opening = KzgEvaluationProof::<E>::new(
            proof.first_msg.f_cm,
            proof.third_msg.f_opening.1,
            x,
            proof.third_msg.f_opening.0,
        );

        let m_opening = KzgEvaluationProof::<E>::new(
            proof.first_msg.m_cm,
            proof.third_msg.m_opening.1,
            x,
            proof.third_msg.m_opening.0,
        );

        let t_opening = KzgEvaluationProof::<E>::new(
            vk.table_cm,
            proof.third_msg.t_opening.1,
            x,
            proof.third_msg.t_opening.0,
        );

        let phi_opening = KzgEvaluationProof::<E>::new(
            proof.second_msg.phi_cm,
            proof.third_msg.phi_opening.1,
            x,
            proof.third_msg.phi_opening.0,
        );

        let root_of_unity = domain.element(1);
        let phi_shifted_opening = KzgEvaluationProof::<E>::new(
            proof.second_msg.phi_cm,
            proof.third_msg.phi_shifted_opening.1,
            x * root_of_unity,
            proof.third_msg.phi_shifted_opening.0,
        );

        let (lhs, rhs) = batch_pairings(
            &[
                f_opening,
                m_opening,
                t_opening,
                phi_opening,
                phi_shifted_opening,
            ],
            r,
        );

        let pairing_check = E::product_of_pairings(&[
            (lhs.into(), vk.g2_gen.clone()),
            (rhs.neg().into(), vk.x.clone()),
        ]);

        if pairing_check != E::Fqk::one() {
            return Err(Error::MVPairingFailed);
        }

        Ok(())
    }
}
