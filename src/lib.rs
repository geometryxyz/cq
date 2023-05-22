mod batcher;
pub mod data_structures;
pub mod error;
pub mod indexer;
pub mod kzg;
pub mod prover;
pub mod rng;
pub mod table;
pub mod tools;
pub mod transcript;
pub mod utils;
pub mod verifier;

pub const PROTOCOL_NAME: &[u8] = b"CQ-1.0";

#[cfg(test)]
mod roundtrip_test {
    use ark_bn254::Bn254;
    use ark_ec::PairingEngine;
    use ark_std::{
        rand::{rngs::StdRng, Rng, RngCore},
        test_rng, UniformRand,
    };
    use rand_chacha::ChaChaRng;
    use sha3::Keccak256;

    use crate::{
        data_structures::{ProvingKey, Statement, Witness},
        indexer::{CommonPreprocessedInput, Index},
        kzg::Kzg,
        prover::Prover,
        rng::SimpleHashFiatShamirRng,
        table::Table,
        utils::unsafe_setup_from_rng,
        verifier::{Verifier, VerifierKey},
    };

    type FS = SimpleHashFiatShamirRng<Keccak256, ChaChaRng>;
    type PrepareResult<E> = (
        Table<<E as PairingEngine>::Fr>,
        Index<E>,
        Statement<E>,
        CommonPreprocessedInput<E>,
        ProvingKey<E>,
        VerifierKey<E>,
        Witness<<E as PairingEngine>::Fr>,
    );

    fn prepare<E: PairingEngine, R: RngCore>(
        n: usize,
        subvector_indices: &[usize],
        rng: &mut R,
    ) -> PrepareResult<E> {
        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<E, R>(n - 1, n, rng);
        let pk = ProvingKey::<E> { srs_g1 };

        let table_values: Vec<_> = (0..n).map(|_| E::Fr::rand(rng)).collect();
        let table = Table::new(&table_values).unwrap();

        let index = Index::<E>::gen(&pk.srs_g1, &srs_g2, &table);

        let witness_values: Vec<_> = subvector_indices.iter().map(|&i| table_values[i]).collect();
        let witness = Witness::<E::Fr>::new(&witness_values).unwrap();

        let statement = Statement::<E> {
            f: Kzg::<E>::commit_g1(&pk.srs_g1, &witness.f).into(),
        };

        let vk = VerifierKey::<E>::new(&srs_g2, table.size, witness.size);
        let common = Index::<E>::compute_common(&srs_g2, &table);

        (table, index, statement, common, pk, vk, witness)
    }

    #[test]
    fn test_roundtrip() {
        let two: usize = 2;
        let n = two.pow(6);

        let mut rng = test_rng();

        let witness_size = two.pow(3);
        let subvector_indices: Vec<usize> =
            (0..witness_size).map(|_| rng.gen_range(0..n - 1)).collect();

        let (table, index, statement, common, pk, vk, witness) =
            prepare::<Bn254, StdRng>(n, &subvector_indices, &mut rng);

        let proof = Prover::<Bn254, FS>::prove(&pk, &index, &table, &witness, &statement).unwrap();
        let res = Verifier::<Bn254, FS>::verify(&vk, &common, &statement, &proof);
        assert!(res.is_ok());
    }
}
