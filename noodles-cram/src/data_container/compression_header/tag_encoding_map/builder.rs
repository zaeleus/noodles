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
        for (tag, value) in record.tags().iter() {
            let key = Key::new(tag, value.ty());
            self.keys.insert(key);
        }
    }

    pub fn build(self) -> TagEncodingMap {
        let mut map = HashMap::new();

        for key in self.keys {
            let id = block::ContentId::from(key);

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
    use noodles_sam::alignment::{
        record::data::field::{Tag, Type},
        record_buf::data::field::Value,
    };

    use super::*;

    #[test]
    fn test_build() {
        let mut builder = Builder::default();

        let mut record = Record::default();
        record.tags.insert(Tag::ALIGNMENT_HIT_COUNT, Value::Int8(1));
        builder.update(&record);

        let mut record = Record::default();
        record.tags.insert(Tag::ALIGNMENT_HIT_COUNT, Value::Int8(1));
        builder.update(&record);

        let mut record = Record::default();
        record.tags.insert(Tag::COMMENT, Value::from("noodles"));
        builder.update(&record);

        let actual = builder.build();

        let nh = block::ContentId::from(Key::new(Tag::ALIGNMENT_HIT_COUNT, Type::Int8));
        let co = block::ContentId::from(Key::new(Tag::COMMENT, Type::String));

        let expected = [
            (
                nh,
                Encoding::new(ByteArray::ByteArrayLen(
                    Encoding::new(Integer::External(nh)),
                    Encoding::new(Byte::External(nh)),
                )),
            ),
            (
                co,
                Encoding::new(ByteArray::ByteArrayLen(
                    Encoding::new(Integer::External(co)),
                    Encoding::new(Byte::External(co)),
                )),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(*actual, expected);
    }
}
