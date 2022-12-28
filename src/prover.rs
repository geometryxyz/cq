use std::{marker::PhantomData, collections::BTreeMap};

use ark_ec::PairingEngine;

use crate::{indexer::Index, data_structures::Witness};

pub struct Prover<E: PairingEngine> {
    _e: PhantomData<E>
}

pub struct State<'a, E: PairingEngine> {
    index: &'a Index<E>,
    witness: &'a Witness<E::Fr>,
}

impl<E: PairingEngine> Prover<E> {
    pub fn prove(

    ) {

    }

    fn round_1<'a>(state: &'a State<E>) {
        let index_repetitions_mapping = BTreeMap::<usize, usize>::default();

        for fi in &state.witness.f_evals {
            // let mut num_of_repetitions = index_repetitions_mapping.entry(fi).or_insert(0);
        }
    }
}