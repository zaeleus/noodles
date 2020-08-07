use std::{collections::HashMap, ops::Deref};

use crate::{record, Record};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagIdsDictionary(Vec<Vec<record::tag::Key>>);

impl TagIdsDictionary {
    pub fn from_records(records: &[Record]) -> Self {
        let mut keys_indices = HashMap::new();

        for record in records {
            let keys: Vec<_> = record.tags.iter().map(|tag| tag.key()).collect();
            let next_index = keys_indices.len();
            keys_indices.entry(keys).or_insert(next_index);
        }

        let mut lines: Vec<_> = keys_indices.into_iter().map(|entry| entry).collect();
        lines.sort_by_key(|(_, index)| *index);
        let dictionary = lines.into_iter().map(|(keys, _)| keys).collect();

        Self(dictionary)
    }
}

impl Deref for TagIdsDictionary {
    type Target = [Vec<record::tag::Key>];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<Vec<Vec<record::tag::Key>>> for TagIdsDictionary {
    fn from(dictionary: Vec<Vec<record::tag::Key>>) -> Self {
        Self(dictionary)
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam::record::data::field::{value::Type, Value};

    use crate::record::{tag::Key, Tag};

    use super::*;

    #[test]
    fn test_from_records() {
        let mut records = Vec::with_capacity(3);

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(1)));
        records.push(record);

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(2)));
        records.push(record);

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(1)));
        record.add_tag(Tag::new(
            Key::new([b'C', b'O'], Type::String),
            Value::String(String::from("noodles")),
        ));
        records.push(record);

        let dictionary = TagIdsDictionary::from_records(&records);
        assert_eq!(
            *dictionary,
            [
                vec![Key::new([b'N', b'H'], Type::Int8)],
                vec![
                    Key::new([b'N', b'H'], Type::Int8),
                    Key::new([b'C', b'O'], Type::String)
                ]
            ]
        )
    }
}
