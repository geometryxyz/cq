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
        prover::{Prover, ProverFirstMessage, ProverSecondMessage, ProverThirdMessage},
        rng::SimpleHashFiatShamirRng,
        table::Table,
        utils::{to_field, unsafe_setup_from_rng},
        verifier::{Verifier, VerifierKey},
    };

    type FS = SimpleHashFiatShamirRng<Keccak256, ChaChaRng>;

    #[test]
    fn test_full_protocol() {
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

        let proof = Prover::<Bn254, FS>::prove(&pk, &index, &table, &witness, &statement).unwrap();

        let vk = VerifierKey::<Bn254>::new(&srs_g2, table.size, witness.size);
        let common = Index::<Bn254>::compute_common(&srs_g2, &table);

        let res = Verifier::<Bn254, FS>::verify(&vk, &common, &statement, &proof);
        assert!(res.is_ok());
    }
}
