use ark_ff::Field;

pub struct Table<F: Field> {
    pub(crate) size: usize,
    pub(crate) values: Vec<F>
}