use std::collections::{HashMap, HashSet};

use super::TagEncodingMap;

use crate::{
    data_container::compression_header::{preservation_map::tag_ids_dictionary::Key, Encoding},
    Record,
};

#[derive(Debug, Default)]
pub struct Builder {
    keys: HashSet<Key>,
}

impl Builder {
    pub fn update(&mut self, record: &Record) {
        for field in record.tags().values() {
            self.keys.insert(field.into());
        }
    }

    pub fn build(self) -> TagEncodingMap {
        let mut map = HashMap::new();

        for key in self.keys {
            let id = key.id();

            let len_encoding = Encoding::External(id);
            let value_encoding = Encoding::External(id);
            let encoding = Encoding::ByteArrayLen(Box::new(len_encoding), Box::new(value_encoding));

            map.insert(id, encoding);
        }

        TagEncodingMap::from(map)
    }
}

#[cfg(test)]
mod tests {
    use noodles_sam::record::data::{
        field::{value::Type, Tag, Value},
        Field,
    };

    use super::*;

    #[test]
    fn test_build() {
        let mut builder = Builder::default();

        let mut record = Record::default();
        record
            .tags
            .insert(Field::new(Tag::AlignmentHitCount, Value::Int8(1)));
        builder.update(&record);

        let mut record = Record::default();
        record
            .tags
            .insert(Field::new(Tag::AlignmentHitCount, Value::Int8(1)));
        builder.update(&record);

        let mut record = Record::default();
        record.tags.insert(Field::new(
            Tag::Comment,
            Value::String(String::from("noodles")),
        ));
        builder.update(&record);

        let actual = builder.build();

        let nh = Key::new(Tag::AlignmentHitCount, Type::Int8);
        let co = Key::new(Tag::Comment, Type::String);

        let expected = [
            (
                nh.id(),
                Encoding::ByteArrayLen(
                    Box::new(Encoding::External(nh.id())),
                    Box::new(Encoding::External(nh.id())),
                ),
            ),
            (
                co.id(),
                Encoding::ByteArrayLen(
                    Box::new(Encoding::External(co.id())),
                    Box::new(Encoding::External(co.id())),
                ),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(*actual, expected);
    }
}
