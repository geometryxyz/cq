pub mod data_structures;
pub mod prover;
pub mod transcript;
pub mod verifier;

pub const MV_PROTOCOL_NAME: &[u8] = b"MV-1.0";
use prover::*;
use verifier::*;

#[cfg(test)]
mod tests {
    use crate::tools::compute_lagrange_basis_commitments;
    use crate::{
        data_structures::Witness,
        kzg::Kzg,
        mv::{data_structures::Statement, prover::Prover, verifier::Verifier},
        rng::SimpleHashFiatShamirRng,
        table::Table,
        utils::{to_field, unsafe_setup_from_rng},
    };
    use ark_ec::{msm::VariableBaseMSM, PairingEngine};
    use ark_ff::PrimeField;
    use ark_poly::{GeneralEvaluationDomain, EvaluationDomain};
    use ark_poly::{univariate::DensePolynomial, UVPolynomial};
    use ark_std::{rand::RngCore, test_rng};
    use rand_chacha::ChaChaRng;
    use sha3::Keccak256;

    use super::{prover::ProvingKey, verifier::VerifierKey};

    fn prepare<E: PairingEngine, R: RngCore>(
        table_values: &Vec<E::Fr>,
        rng: &mut R,
    ) -> (ProvingKey<E>, VerifierKey<E>) {
        let n = table_values.len();
        let (srs_g1, srs_g2) = unsafe_setup_from_rng::<E, R>(2*n - 1, 1, rng);

        let pk = ProvingKey::<E> { srs_g1 };

        let domain = GeneralEvaluationDomain::<E::Fr>::new(n).unwrap();
        let table_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&table_values));
        let table_cm = Kzg::<E>::commit_g1(&pk.srs_g1, &table_poly).into();

        // let table_values = table_values
        //     .iter()
        //     .map(|ti| ti.into_repr())
        //     .collect::<Vec<_>>();
        // let table_cm = VariableBaseMSM::multi_scalar_mul(&li_commitments, &table_values);

        let vk = VerifierKey::<E>::new(&srs_g2, n, table_cm.into());

        (pk, vk)
    }

    #[test]
    fn test_verifier() {
        use ark_bn254::{Bn254, Fq, Fr, G1Affine, G2Affine};
        type FS = SimpleHashFiatShamirRng<Keccak256, ChaChaRng>;

        let n = 8;
        let domain = GeneralEvaluationDomain::<Fr>::new(n).unwrap();
        let table_values = vec![1, 5, 10, 15, 20, 25, 30, 100];
        let table = Table::new(&to_field(&table_values)).unwrap();
        let table_poly = DensePolynomial::from_coefficients_slice(&domain.ifft(&table.values));

        let witness_values = vec![5, 5, 20, 1, 100, 100, 15, 100];
        let witness = Witness::<Fr>::new(&to_field(&witness_values)).unwrap();

        let mut rng = test_rng();

        let (pk, vk) = prepare(&table.values, &mut rng);

        let statement = Statement::<Bn254> {
            f_cm: Kzg::<Bn254>::commit_g1(&pk.srs_g1, &witness.f).into(),
        };

        let proof = Prover::<Bn254, FS>::prove(&pk, &vk, &table, &witness, &statement, &table_poly)
            .unwrap();
        let res = Verifier::<Bn254, FS>::verify(&vk, &statement, &proof).unwrap();

        // assert!(res.is_ok());
    }
}
