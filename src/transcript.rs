use ark_ec::PairingEngine;
use ark_ff::{ToBytes, to_bytes};

use crate::{rng::FiatShamirRng, prover::ProverFirstMessage};

pub struct TranscriptOracle<FS: FiatShamirRng> {
    fs_rng: FS
}

impl<FS: FiatShamirRng> TranscriptOracle<FS> {
    pub fn initialize<'a, T: 'a + ToBytes>(initial_input: &'a T) -> Self {
        let fs_rng = FS::initialize(&to_bytes![initial_input].unwrap());
        Self {
            fs_rng
        }
    }

    pub fn stream_first_message<E: PairingEngine>(&mut self, msg: &ProverFirstMessage<E>) {
        self.fs_rng.absorb(&to_bytes![msg].unwrap());
    }
}