use std::{
    collections::{HashMap, HashSet},
    ops::Deref,
};

use crate::{
    num::{write_itf8, Itf8},
    Record,
};

use super::{encoding, Encoding};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TagEncodingMap(HashMap<Itf8, Encoding>);

impl TagEncodingMap {
    pub fn from_records(records: &[Record]) -> Self {
        let mut keys = HashSet::new();

        for record in records {
            for tag in &record.tags {
                keys.insert(tag.key());
            }
        }

        let mut map = HashMap::new();

        for key in keys {
            let id = key.id();

            let mut args = Vec::new();
            write_itf8(&mut args, id).unwrap();

            // TODO: Select encoding depending on the type of data.
            let encoding = Encoding::new(encoding::Kind::External, args);

            map.insert(id, encoding);
        }

        Self::from(map)
    }
}

impl Deref for TagEncodingMap {
    type Target = HashMap<Itf8, Encoding>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<HashMap<Itf8, Encoding>> for TagEncodingMap {
    fn from(map: HashMap<Itf8, Encoding>) -> Self {
        Self(map)
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

        let nh = Key::new([b'N', b'H'], Type::Int8);
        let co = Key::new([b'C', b'O'], Type::String);

        let mut record = Record::default();
        record.add_tag(Tag::new(nh, Value::Int8(1)));
        records.push(record);

        let mut record = Record::default();
        record.add_tag(Tag::new(nh, Value::Int8(1)));
        records.push(record);

        let mut record = Record::default();
        record.add_tag(Tag::new(co, Value::String(String::from("noodles"))));
        records.push(record);

        let actual = TagEncodingMap::from_records(&records);

        let expected = vec![
            (
                nh.id(),
                Encoding::new(encoding::Kind::External, vec![0xe0, 0x4e, 0x48, 0x63]),
            ),
            (
                co.id(),
                Encoding::new(encoding::Kind::External, vec![0xe0, 0x43, 0x4f, 0x5a]),
            ),
        ]
        .into_iter()
        .collect();

        assert_eq!(*actual, expected);
    }
}
