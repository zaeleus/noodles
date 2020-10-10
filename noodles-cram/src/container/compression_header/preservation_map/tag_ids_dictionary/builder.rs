use std::collections::HashMap;

use crate::{record, Record};

use super::TagIdsDictionary;

#[derive(Debug, Default)]
pub struct Builder {
    keys_indices: HashMap<Vec<record::tag::Key>, usize>,
}

impl Builder {
    pub fn update(&mut self, record: &Record) {
        let keys: Vec<_> = record.tags().iter().map(|tag| tag.key()).collect();
        let next_index = self.keys_indices.len();
        self.keys_indices.entry(keys).or_insert(next_index);
    }

    pub fn build(self) -> TagIdsDictionary {
        let mut lines: Vec<_> = self.keys_indices.into_iter().collect();
        lines.sort_by_key(|(_, index)| *index);
        let dictionary: Vec<_> = lines.into_iter().map(|(keys, _)| keys).collect();
        TagIdsDictionary::from(dictionary)
    }
}

#[cfg(test)]
mod tests {
    use noodles_bam::record::data::field::{value::Type, Value};

    use crate::record::{tag::Key, Tag};

    use super::*;

    #[test]
    fn test_from_records() {
        let mut builder = Builder::default();

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(1)));
        builder.update(&record);

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(2)));
        builder.update(&record);

        let mut record = Record::default();
        record.add_tag(Tag::new(Key::new([b'N', b'H'], Type::Int8), Value::Int8(1)));
        record.add_tag(Tag::new(
            Key::new([b'C', b'O'], Type::String),
            Value::String(String::from("noodles")),
        ));
        builder.update(&record);

        let dictionary = builder.build();

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
