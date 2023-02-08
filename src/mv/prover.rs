use std::marker::PhantomData;

use ark_ec::PairingEngine;
use ark_ff::{batch_inversion, ToBytes, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use ark_std::{One, Zero};

use crate::{
    data_structures::Witness, error::Error, kzg::Kzg, mv::MV_PROTOCOL_NAME, rng::FiatShamirRng,
    table::Table,
};

use super::{data_structures::Statement, transcript::TranscriptOracle, verifier::VerifierKey};

pub struct ProvingKey<E: PairingEngine> {
    pub(crate) srs_g1: Vec<E::G1Affine>,
}

pub struct ProverFirstMessage<E: PairingEngine> {
    /// Commitment to the input polynomial
    pub(super) f_cm: E::G1Affine,
    /// Commitment to the multiplicities polynomial
    pub(super) m_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for ProverFirstMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.f_cm.write(&mut w)?;
        self.m_cm.write(&mut w)
    }
}

pub struct ProverSecondMessage<E: PairingEngine> {
    /// Commitment to grand sum polynomial
    pub(super) phi_cm: E::G1Affine,
    /// Commitment to quotient polynomial
    pub(super) q_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for ProverSecondMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.phi_cm.write(&mut w)
    }
}

pub struct ProverThirdMessage<E: PairingEngine> {
    /// Opening of f polynomial commitment
    pub(super) f_opening: (E::Fr, E::G1Affine),
    /// Opening of m polynomial commitment
    pub(super) m_opening: (E::Fr, E::G1Affine),
    /// Opening of table
    pub(super) t_opening: (E::Fr, E::G1Affine),
    /// Opening of grand sum polynomial commitment at x
    pub(super) phi_opening: (E::Fr, E::G1Affine),
    /// Opening of grand sum polynomial commitment at omega * x
    pub(super) phi_shifted_opening: (E::Fr, E::G1Affine),
    /// Opening of quotient polynomial
    pub(super) q_opening: (E::Fr, E::G1Affine),
}

impl<E: PairingEngine> ToBytes for ProverThirdMessage<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.m_opening.0.write(&mut w)?;
        self.m_opening.1.write(&mut w)?;
        self.phi_opening.0.write(&mut w)?;
        self.phi_opening.1.write(&mut w)?;
        self.t_opening.0.write(&mut w)?;
        self.t_opening.1.write(&mut w)?;
        self.phi_shifted_opening.0.write(&mut w)?;
        self.phi_shifted_opening.1.write(&mut w)?;
        self.q_opening.0.write(&mut w)?;
        self.q_opening.1.write(&mut w)
    }
}

pub struct Prover<E: PairingEngine, FS: FiatShamirRng>(PhantomData<(E, FS)>);

pub struct State<'a, E: PairingEngine> {
    pk: &'a ProvingKey<E>,
    table: &'a Table<E::Fr>,
    witness: &'a Witness<E::Fr>,

    f_poly: Option<DensePolynomial<E::Fr>>,
    table_poly: &'a DensePolynomial<E::Fr>,

    /// Multiplicities of table elements in input
    m_values: Option<Vec<E::Fr>>,
    m_poly: Option<DensePolynomial<E::Fr>>,

    /// grand sum values
    phi_values: Option<Vec<E::Fr>>,
    phi_poly: Option<DensePolynomial<E::Fr>>,

    // quotient polynomial
    q_poly: Option<DensePolynomial<E::Fr>>,
}

impl<'a, E: PairingEngine> State<'a, E> {
    pub fn new(
        pk: &'a ProvingKey<E>,
        table: &'a Table<E::Fr>,
        witness: &'a Witness<E::Fr>,
        table_poly: &'a DensePolynomial<E::Fr>,
    ) -> Self {
        Self {
            pk,
            table,
            witness,
            f_poly: None,
            table_poly,
            m_values: None,
            m_poly: None,
            phi_values: None,
            phi_poly: None,
            q_poly: None,
        }
    }
}

pub struct Proof<E: PairingEngine> {
    pub(crate) first_msg: ProverFirstMessage<E>,
    pub(crate) second_msg: ProverSecondMessage<E>,
    pub(crate) third_msg: ProverThirdMessage<E>,
}

impl<E: PairingEngine, FS: FiatShamirRng> Prover<E, FS> {
    pub fn prove<'a>(
        pk: &'a ProvingKey<E>,
        vk: &'a VerifierKey<E>,
        table: &'a Table<E::Fr>,
        witness: &'a Witness<E::Fr>,
        statement: &Statement<E>,
        table_poly: &'a DensePolynomial<E::Fr>,
    ) -> Result<Proof<E>, Error> {
        let mut state = State::new(pk, table, witness, table_poly);
        let mut transcript = TranscriptOracle::<FS>::initialize(&MV_PROTOCOL_NAME);

        transcript.stream_public_input(vk, statement);

        let first_msg = Self::round_1(&mut state)?;
        transcript.stream_first_message(&first_msg);

        let beta: E::Fr = transcript.squeeze_challenge();

        let second_msg = Self::round_2(&mut state, beta)?;
        transcript.stream_second_message(&second_msg);

        let x: E::Fr = transcript.squeeze_challenge();

        let third_msg = Self::round_3(&mut state, x)?;

        Ok(Proof {
            first_msg,
            second_msg,
            third_msg,
        })
    }

    pub fn round_1(state: &mut State<E>) -> Result<ProverFirstMessage<E>, Error> {
        let table_n = state.table.size;
        let domain = GeneralEvaluationDomain::<E::Fr>::new(table_n).unwrap();

        let mut m_values = vec![E::Fr::zero(); table_n];

        // We expect that witness is padded to table_n size
        assert_eq!(state.witness.f_evals.len(), table_n);

        for fi in &state.witness.f_evals {
            // Find index of f_i in table
            let index = state.table.value_index_mapping.get(fi);
            let index = index.ok_or(Error::ValueNotInTable(format!("{}", fi)))?;

            // Increment the multiplicity of `index`.
            m_values[*index] += E::Fr::one();
        }

        let f_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&state.witness.f_evals));
        let f_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &f_poly).into();
        state.f_poly = Some(f_poly);

        let m_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&m_values));
        let m_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &m_poly).into();
        state.m_values = Some(m_values);
        state.m_poly = Some(m_poly);

        Ok(ProverFirstMessage { f_cm, m_cm })
    }

    pub fn round_2(state: &mut State<E>, beta: E::Fr) -> Result<ProverSecondMessage<E>, Error> {
        let table_n = state.table.size;
        let domain = GeneralEvaluationDomain::<E::Fr>::new(table_n).unwrap();

        let mut phi_values = vec![E::Fr::zero()];

        let mut lhs: Vec<_> = state
            .witness
            .f_evals
            .iter()
            .map(|&f_i| f_i + beta)
            .collect();
        batch_inversion(&mut lhs);

        let mut rhs: Vec<_> = state.table.values.iter().map(|&t_i| t_i + beta).collect();
        batch_inversion(&mut rhs);

        let rhs: Vec<_> = rhs
            .iter()
            .zip(
                state
                    .m_values
                    .as_ref()
                    .expect("m evals not in state")
                    .iter(),
            )
            .map(|(&rh_i, &m_i)| m_i * rh_i)
            .collect();

        for i in 0..table_n {
            // s = (1 / (beta + f_i)) - (m[i] / (beta + t_i))
            let s = lhs[i] - rhs[i];
            phi_values.push(phi_values[i] + s);
        }
        if phi_values[table_n] != E::Fr::zero() {
            return Err(Error::ProverGrandSumFailed);
        }
        phi_values.pop();


        let phi_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&phi_values));
        state.phi_values = Some(phi_values);
        let phi_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &phi_poly).into();

        //Compute quotient
        let m_poly = state.m_poly.as_ref().expect("m poly not in the state");
        let f_poly = state.f_poly.as_ref().expect("f poly not in the state");

        let table_poly = state.table_poly;

        let domain_2n = GeneralEvaluationDomain::<E::Fr>::new(table_n * 2).unwrap();

        let extended_m_values = domain_2n.coset_fft(m_poly);
        let extended_f_values = domain_2n.coset_fft(f_poly);

        // TODO: put this in table
        let extended_table_values = domain_2n.coset_fft(table_poly);
        let extended_phi_values = domain_2n.coset_fft(&phi_poly);

        let next_rotation_in_2n = |i: usize| extended_phi_values[(i + 2) % domain_2n.size()];

        let mut extended_zh_values = compute_vanishing_poly_over_coset(domain_2n, table_n as u64);
        batch_inversion(&mut extended_zh_values);

        // phi(Xg) - phi(X) = 1/(beta + f) - m/(beta + t)
        // (beta + t) * (beta + f) * (phi(Xg) - phi(X)) - ((beta + t) - m*(beta + f)) =  Q(X) * ZH
        let mut q_values = vec![E::Fr::zero(); domain_2n.size()];
        for i in 0..domain_2n.size() {
            let lhs = {
                (beta + extended_table_values[i])
                    * (beta + extended_f_values[i])
                    * (next_rotation_in_2n(i) - extended_phi_values[i])
            };
            let rhs = {
                beta + extended_table_values[i] - extended_m_values[i] * (beta + extended_f_values[i])
            };
            q_values[i] = (lhs - rhs) * extended_zh_values[i];
        }

        let q_poly = DensePolynomial::from_coefficients_slice(&domain_2n.coset_ifft(&q_values));
        let q_cm: E::G1Affine = Kzg::<E>::commit_g1(&state.pk.srs_g1, &q_poly).into();

        state.phi_poly = Some(phi_poly);
        state.q_poly = Some(q_poly);

        Ok(ProverSecondMessage { phi_cm, q_cm })
    }

    pub fn round_3(state: &mut State<E>, x: E::Fr) -> Result<ProverThirdMessage<E>, Error> {
        let f_poly = state.f_poly.as_ref().expect("f poly not in the state");
        let m_poly = state.m_poly.as_ref().expect("m poly not in the state");
        let phi_poly = state.phi_poly.as_ref().expect("phi poly not in the state");
        let table_poly = state.table_poly;
        let q_poly = state.q_poly.as_ref().expect("q poly poly not in the state");

        let table_n = state.table.size;
        let domain = GeneralEvaluationDomain::<E::Fr>::new(table_n).unwrap();

        let f_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, f_poly, x);
        let m_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, m_poly, x);

        let t_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, &table_poly, x);

        let phi_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, phi_poly, x);

        let root_of_unity = domain.element(1);
        let phi_shifted_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, phi_poly, x * root_of_unity);

        let q_opening = Kzg::<E>::open_g1(&state.pk.srs_g1, q_poly, x);

        Ok(ProverThirdMessage {
            f_opening,
            m_opening,
            t_opening,
            phi_opening,
            phi_shifted_opening,
            q_opening,
        })
    }
}

pub fn compute_vanishing_poly_over_coset<F, D>(
    domain: D,        // domain to evaluate over
    poly_degree: u64, // degree of the vanishing polynomial
) -> Vec<F>
where
    F: PrimeField,
    D: EvaluationDomain<F>,
{
    assert!(
        (domain.size() as u64) > poly_degree,
        "domain_size = {}, poly_degree = {}",
        domain.size() as u64,
        poly_degree
    );
    let group_gen = domain.element(1);
    let coset_gen = F::multiplicative_generator().pow([poly_degree, 0, 0, 0]);
    let v_h: Vec<_> = (0..domain.size())
        .map(|i| {
            (coset_gen * group_gen.pow([poly_degree * i as u64, 0, 0, 0]))
                - F::one()
        })
        .collect();
    v_h
}
