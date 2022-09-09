use std::collections::{HashMap, HashSet};

use super::TagEncodingMap;

use crate::{
    container::block,
    data_container::compression_header::{
        encoding::codec::{Byte, ByteArray, Integer},
        preservation_map::tag_ids_dictionary::Key,
        Encoding,
    },
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
            let id = block::ContentId::from(key.id());

            let len_encoding = Encoding::new(Integer::External(id));
            let value_encoding = Encoding::new(Byte::External(id));
            let encoding = Encoding::new(ByteArray::ByteArrayLen(len_encoding, value_encoding));

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
                block::ContentId::from(nh.id()),
                Encoding::new(ByteArray::ByteArrayLen(
                    Encoding::new(Integer::External(block::ContentId::from(nh.id()))),
                    Encoding::new(Byte::External(block::ContentId::from(nh.id()))),
                )),
            ),
            (
                block::ContentId::from(co.id()),
                Encoding::new(ByteArray::ByteArrayLen(
                    Encoding::new(Integer::External(block::ContentId::from(co.id()))),
                    Encoding::new(Byte::External(block::ContentId::from(co.id()))),
                )),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(*actual, expected);
    }
}
