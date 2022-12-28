use std::collections::BTreeMap;

use ark_ff::FftField;

use crate::error::Error;

#[derive(Debug)]
pub struct Table<F: FftField> {
    pub(crate) size: usize,
    pub(crate) values: Vec<F>,
    pub(crate) value_index_mapping: BTreeMap<F, usize>,
}

impl<F: FftField> Table<F> {
    pub fn new(values: &Vec<F>) -> Result<Self, Error> {
        if !values.len().is_power_of_two() {
            return Err(Error::TableSizeNotPow2(values.len()));
        }
        let mut value_index_mapping = BTreeMap::<F, usize>::default();
        for (i, &ti) in values.iter().enumerate() {
            let prev = value_index_mapping.insert(ti, i);
            if prev.is_some() {
                return Err(Error::DuplicateValueInTable(format!("{}", ti)));
            }
        }
        Ok(Self {
            size: values.len(),
            values: values.clone(),
            value_index_mapping,
        })
    }
}

#[cfg(test)]
pub mod table_tests {
    use crate::error::Error;
    use ark_bn254::Fr;
    use ark_ff::UniformRand;
    use ark_std::test_rng;

    use super::Table;

    #[test]
    fn test_correct_table() {
        let n = 32;
        let mut rng = test_rng();

        let table_values: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let _ = Table::new(&table_values);
    }

    #[test]
    fn test_not_pow_2() {
        let n = 31;
        let mut rng = test_rng();

        let table_values: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        let res = Table::new(&table_values);
        assert_eq!(res.unwrap_err(), Error::TableSizeNotPow2(n));
    }

    #[test]
    fn test_dup_value() {
        let n = 32;
        let mut rng = test_rng();

        let mut table_values: Vec<_> = (0..n).map(|_| Fr::rand(&mut rng)).collect();
        table_values[5] = table_values[10];

        let res = Table::new(&table_values);
        assert_eq!(
            res.unwrap_err(),
            Error::DuplicateValueInTable(format!("{}", table_values[5]))
        );
    }
}
