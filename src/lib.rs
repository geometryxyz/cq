pub mod data_structures;
pub mod error;
pub mod indexer;
pub mod kzg;
pub mod prover;
pub mod table;
pub mod tools;
pub mod utils;
pub mod rng;
pub mod transcript;
pub mod verifier;

pub const PROTOCOL_NAME: &'static [u8] = b"CQ-1.0";

#[cfg(test)]
mod roundtrip_test {
    use std::ops::Neg;

    use ark_bn254::{Bn254, Fr, G2Affine, Fq12, G1Affine};
    use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
    use ark_ff::{UniformRand, One, Field};
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use ark_std::{rand::rngs::StdRng, test_rng};
    use rand_chacha::ChaChaRng;
    use sha3::Keccak256;

    use crate::{
        data_structures::{ProvingKey, Witness, Statement},
        indexer::Index,
        table::Table,
        utils::{to_field, unsafe_setup_from_rng}, kzg::Kzg, rng::SimpleHashFiatShamirRng,
        prover::{Prover, ProverFirstMessage, ProverSecondMessage, ProverThirdMessage}, verifier::{VerifierKey, Verifier}
    };

    type FS = SimpleHashFiatShamirRng<Keccak256, ChaChaRng>;

    #[test]
    fn test_full_protocol() {
        let n = 8;
        let mut rng = test_rng();

        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<Bn254, StdRng>(n - 1, n, &mut rng);
        let pk = ProvingKey { srs_g1, srs_g2 };

        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 35];
        let table = Table::new(&to_field(&table_values)).unwrap();

        let index = Index::<Bn254>::gen(&pk.srs_g1, &pk.srs_g2, &table);

        let witness_values = vec![5, 15, 20, 35];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let statement = Statement::<Bn254> {
            f: Kzg::<Bn254>::commit_g1(&pk.srs_g1, &witness.f).into()
        };

        let proof = Prover::<Bn254, FS>::prove(&pk, &index, &table, &witness, &statement).unwrap();

        let vk = VerifierKey::<Bn254>::new(&pk.srs_g2, table.size, witness.size);
        let common = Index::<Bn254>::compute_common(&pk.srs_g2, &table);

        let res = Verifier::<Bn254, FS>::verify(&vk, &common, &statement, &proof);
        assert!(res.is_ok());
    }
}