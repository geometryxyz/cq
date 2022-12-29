use ark_ec::PairingEngine;
use ark_ff::{to_bytes, Field, ToBytes};

use crate::{
    data_structures::Statement,
    indexer::{CommonPreprocessedInput, Index},
    prover::{ProverFirstMessage, ProverSecondMessage, ProverThirdMessage},
    rng::FiatShamirRng,
    table::Table,
};

pub struct TranscriptOracle<FS: FiatShamirRng> {
    fs_rng: FS,
}

impl<FS: FiatShamirRng> TranscriptOracle<FS> {
    pub fn initialize<'a, T: 'a + ToBytes>(initial_input: &'a T) -> Self {
        let fs_rng = FS::initialize(&to_bytes![initial_input].unwrap());
        Self { fs_rng }
    }

    pub fn squeeze_challenge<F: Field>(&mut self) -> F {
        F::rand(&mut self.fs_rng)
    }

    pub fn stream_public_input<E: PairingEngine>(
        &mut self,
        common: &CommonPreprocessedInput<E>,
        statement: &Statement<E>,
    ) {
        self.fs_rng.absorb(&to_bytes![common, statement].unwrap());
    }

    pub fn stream_first_message<E: PairingEngine>(&mut self, msg: &ProverFirstMessage<E>) {
        self.fs_rng.absorb(&to_bytes![msg].unwrap());
    }

    pub fn stream_second_message<E: PairingEngine>(&mut self, msg: &ProverSecondMessage<E>) {
        self.fs_rng.absorb(&to_bytes![msg].unwrap());
    }

    pub fn stream_third_message<E: PairingEngine>(&mut self, msg: &ProverThirdMessage<E>) {
        self.fs_rng.absorb(&to_bytes![msg].unwrap());
    }
}
