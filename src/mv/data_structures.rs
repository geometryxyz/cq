use ark_ec::PairingEngine;
use ark_ff::ToBytes;

pub struct Statement<E: PairingEngine> {
    /// Commitment to witness polynomial
    pub(crate) f_cm: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for Statement<E> {
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        self.f_cm.write(&mut w)
    }
}
