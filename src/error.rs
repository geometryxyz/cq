#[derive(Debug, PartialEq)]
pub enum Error {
    TableSizeNotPow2(usize),
    DuplicateValueInTable(String),
    ValueNotInTable(String)
}