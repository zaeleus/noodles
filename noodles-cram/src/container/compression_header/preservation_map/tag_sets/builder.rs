use std::collections::HashMap;

use super::{Key, TagSets};
use crate::io::writer::Record;

#[derive(Debug, Default)]
pub struct Builder {
    keys_indices: HashMap<Vec<Key>, usize>,
}

impl Builder {
    pub fn update(&mut self, record: &Record) {
        let keys: Vec<_> = record
            .data
            .iter()
            .map(|(tag, value)| Key::new(*tag, value.ty()))
            .collect();

        let next_index = self.keys_indices.len();

        self.keys_indices.entry(keys).or_insert(next_index);
    }

    pub(crate) fn build(self) -> TagSets {
        let mut lines: Vec<_> = self.keys_indices.into_iter().collect();
        lines.sort_by_key(|(_, index)| *index);
        let sets: Vec<_> = lines.into_iter().map(|(keys, _)| keys).collect();
        TagSets::from(sets)
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::alignment::{
        record::data::field::{Tag, Type},
        record_buf::data::field::Value,
    };

    use super::*;

    #[test]
    fn test_from_records() {
        let mut builder = Builder::default();

        let mut record = Record::default();
        record.data.push((Tag::ALIGNMENT_HIT_COUNT, Value::Int8(1)));
        builder.update(&record);

        let mut record = Record::default();
        record.data.push((Tag::ALIGNMENT_HIT_COUNT, Value::Int8(2)));
        builder.update(&record);

        let mut record = Record::default();
        record.data.push((Tag::ALIGNMENT_HIT_COUNT, Value::Int8(1)));
        record.data.push((Tag::COMMENT, Value::from("noodles")));
        builder.update(&record);

        let dictionary = builder.build();

        assert_eq!(
            *dictionary,
            [
                vec![Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8)],
                vec![
                    Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8),
                    Key::new(Tag::COMMENT, Type::String)
                ]
            ]
        )
    }
}
