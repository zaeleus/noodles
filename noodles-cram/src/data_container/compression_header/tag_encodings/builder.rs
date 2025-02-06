use std::collections::HashSet;

use super::TagEncodings;
use crate::{
    container::block,
    data_container::compression_header::{
        encoding::codec::{Byte, ByteArray, Integer},
        preservation_map::tag_sets::Key,
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

    pub fn build(self) -> TagEncodings {
        let mut encodings = TagEncodings::new();

        for key in self.keys {
            let block_content_id = block::ContentId::from(key);

            let len_encoding = Encoding::new(Integer::External { block_content_id });
            let value_encoding = Encoding::new(Byte::External { block_content_id });
            let encoding = Encoding::new(ByteArray::ByteArrayLen {
                len_encoding,
                value_encoding,
            });

            encodings.insert(block_content_id, encoding);
        }

        encodings
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
                Encoding::new(ByteArray::ByteArrayLen {
                    len_encoding: Encoding::new(Integer::External {
                        block_content_id: nh,
                    }),
                    value_encoding: Encoding::new(Byte::External {
                        block_content_id: nh,
                    }),
                }),
            ),
            (
                co,
                Encoding::new(ByteArray::ByteArrayLen {
                    len_encoding: Encoding::new(Integer::External {
                        block_content_id: co,
                    }),
                    value_encoding: Encoding::new(Byte::External {
                        block_content_id: co,
                    }),
                }),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(actual, expected);
    }
}
