#[derive(Debug, PartialEq)]
pub enum Error {
    TableSizeNotPow2(usize),
    WitnessSizeNotPow2(usize),
    DuplicateValueInTable(String),
    ValueNotInTable(String),

    BatchedPairingFailed,

    Pairing1Failed,
    Pairing2Failed,
    Pairing3Failed,
    Pairing4Failed,
}
