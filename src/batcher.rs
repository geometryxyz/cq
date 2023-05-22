use std::collections::HashMap;

use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve};
use ark_ff::One;

pub struct PairingBatcher<E: PairingEngine> {
    /// Mapping of all G2 points serialized with correlated G1 points
    g2_to_g1: HashMap<E::G2Projective, E::G1Projective>,
    /// challenge
    challenge: E::Fr,
    /// running challenge
    running_challenge: E::Fr,
}

impl<E: PairingEngine> PairingBatcher<E> {
    pub fn new(challenge: E::Fr) -> Self {
        Self {
            g2_to_g1: HashMap::default(),
            challenge,
            running_challenge: E::Fr::one(),
        }
    }

    /// Adds new pairing equation that needs to be checked
    pub fn add_pairing(&mut self, pairs: &[(E::G1Affine, E::G2Affine)]) {
        let g2_points: Vec<E::G2Projective> = pairs.iter().map(|&(_, g2)| g2.into()).collect();

        let mut is_present: bool = false;
        for g2 in g2_points.iter() {
            if self.g2_to_g1.get(g2).is_some() {
                is_present = true;
                break;
            }
        }

        let g1_points: Vec<E::G1Projective> = if is_present {
            self.running_challenge *= self.challenge;
            pairs
                .iter()
                .map(|&(g1, _)| g1.mul(self.running_challenge))
                .collect()
        } else {
            pairs.iter().map(|pair| pair.0.into()).collect()
        };

        self.update_mapping(&g2_points, &g1_points);
    }

    /// Updates mapping based on pairs that are added
    fn update_mapping(&mut self, g2_points: &[E::G2Projective], g1_points: &[E::G1Projective]) {
        g2_points
            .iter()
            .zip(g1_points.iter())
            .for_each(|(&g2, g1)| {
                self.g2_to_g1
                    .entry(g2)
                    .and_modify(|g1_point| *g1_point += g1)
                    .or_insert(*g1);
            });
    }

    /// Returns output that is ready to be called on MultiMillerLoop
    pub fn finalize(&self) -> Vec<(E::G1Prepared, E::G2Prepared)> {
        self.g2_to_g1
            .iter()
            .map(|(&g2, &g1)| (g1.into_affine().into(), g2.into_affine().into()))
            .collect()
    }
}

#[cfg(test)]
mod test {
    use ark_bn254::{Bn254, Fq12, Fr, G1Affine, G2Affine};
    use ark_ec::{AffineCurve, PairingEngine};
    use ark_ff::Field;
    use ark_std::{test_rng, One, UniformRand};
    use std::ops::Neg;

    use crate::batcher::PairingBatcher;

    #[test]
    fn test() {
        /*
            e(a, b) = e(c, d)
            e(j, b) = e(f, g)
            e(e, d) = e(h, b)
        */

        let mut rng = test_rng();

        let a = Fr::rand(&mut rng);
        let b = Fr::rand(&mut rng);
        let c = Fr::rand(&mut rng);
        let d = a * b * c.inverse().unwrap();
        let f = Fr::rand(&mut rng);
        let j = Fr::rand(&mut rng);
        let g = j * b * f.inverse().unwrap();
        let e = Fr::rand(&mut rng);
        let h = e * d * b.inverse().unwrap();

        let a: G1Affine = (G1Affine::prime_subgroup_generator().mul(a)).into();
        let b: G2Affine = (G2Affine::prime_subgroup_generator().mul(b)).into();
        let c: G1Affine = (G1Affine::prime_subgroup_generator().mul(c)).into();
        let d: G2Affine = (G2Affine::prime_subgroup_generator().mul(d)).into();
        let j: G1Affine = (G1Affine::prime_subgroup_generator().mul(j)).into();
        let f: G1Affine = (G1Affine::prime_subgroup_generator().mul(f)).into();
        let g: G2Affine = (G2Affine::prime_subgroup_generator().mul(g)).into();
        let e: G1Affine = (G1Affine::prime_subgroup_generator().mul(e)).into();
        let h: G1Affine = (G1Affine::prime_subgroup_generator().mul(h)).into();

        // Manual Miller loop
        /*
            e(a, b) = e(c, d)
            e(j, b) = e(f, g)
            e(e, d) = e(h, b)
        */
        {
            let result = {
                Bn254::product_of_pairings(&[
                    (a.into(), b.into()),
                    (c.neg().into(), d.into()),
                    (j.into(), b.into()),
                    (f.neg().into(), g.into()),
                    (e.into(), d.into()),
                    (h.neg().into(), b.into()),
                ])
            };

            assert_eq!(result, Fq12::one());
        }

        {
            // Batched test
            let mut pairing_batcher = PairingBatcher::<Bn254>::new(Fr::rand(&mut rng));

            pairing_batcher.add_pairing(&[(a, b), ((-c), d)]);
            pairing_batcher.add_pairing(&[(j, b), ((-f), g)]);
            pairing_batcher.add_pairing(&[(e, d), ((-h), b)]);

            let tuples = pairing_batcher.finalize();
            let result = Bn254::product_of_pairings(&tuples);
            assert_eq!(result, Fq12::one());
        }
    }
}
