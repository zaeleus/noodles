//! BAM record data.

pub mod field;

use std::{borrow::Borrow, fmt, io, iter};

use noodles_sam::{
    self as sam,
    alignment::record::data::field::{Tag, Value},
};

use self::field::decode_field;

/// BAM record data.
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub(super) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns whether there are any fields.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the value of the given tag.
    pub fn get<K>(&self, tag: &K) -> Option<io::Result<Value<'_>>>
    where
        K: Borrow<[u8; 2]>,
    {
        for result in self.iter() {
            match result {
                Ok((t, value)) => {
                    if &t == tag.borrow() {
                        return Some(Ok(value));
                    }
                }
                Err(e) => return Some(Err(e)),
            };
        }

        None
    }

    /// Returns an iterator over all tag-value pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_ {
        let mut src = self.0;

        iter::from_fn(move || {
            if src.is_empty() {
                None
            } else {
                Some(decode_field(&mut src))
            }
        })
    }
}

impl fmt::Debug for Data<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut formatter = f.debug_map();

        for result in self.iter() {
            let (tag, value) = result.map_err(|_| fmt::Error)?;
            formatter.entry(&tag, &value);
        }

        formatter.finish()
    }
}

impl sam::alignment::record::Data for Data<'_> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }

    fn get(&self, tag: &Tag) -> Option<io::Result<Value<'_>>> {
        self.get(tag)
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Tag, Value<'_>)>> + '_> {
        Box::new(self.iter())
    }
}

impl AsRef<[u8]> for Data<'_> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Data<'a>> for sam::alignment::record_buf::Data {
    type Error = io::Error;

    fn try_from(bam_data: Data<'a>) -> Result<Self, Self::Error> {
        use crate::record::codec::decoder::read_data;

        let mut src = bam_data.0;
        let mut sam_data = Self::default();
        read_data(&mut src, &mut sam_data)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        Ok(sam_data)
    }
}

pub(super) fn get_raw_cigar<'a>(src: &mut &'a [u8]) -> io::Result<Option<&'a [u8]>> {
    use noodles_sam::alignment::record::data::field::Type;

    use self::field::{
        decode_tag, decode_type, decode_value,
        value::array::{decode_raw_array, decode_subtype},
    };

    fn get_array_field<'a>(src: &mut &'a [u8]) -> io::Result<Option<(Tag, &'a [u8])>> {
        let tag = decode_tag(src)?;
        let ty = decode_type(src)?;

        if ty == Type::Array {
            let subtype = decode_subtype(src)?;
            let buf = decode_raw_array(src, subtype)?;
            Ok(Some((tag, buf)))
        } else {
            decode_value(src, ty)?;
            Ok(None)
        }
    }

    while !src.is_empty() {
        if let Some((Tag::CIGAR, buf)) = get_array_field(src)? {
            return Ok(Some(buf));
        }
    }

    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_empty() {
        let data = Data::new(&[]);
        assert!(data.is_empty());

        let data = Data::new(&[b'N', b'H', b'C', 0x01]);
        assert!(!data.is_empty());
    }

    #[test]
    fn test_get() -> io::Result<()> {
        let data = Data::new(&[b'N', b'H', b'C', 0x01]);

        assert!(data.get(&Tag::ALIGNMENT_HIT_COUNT).is_some());
        assert!(data.get(b"NH").is_some());

        assert!(data.get(&Tag::COMMENT).is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let data = Data::new(&[]);
        assert!(data.iter().next().is_none());

        let data = Data::new(&[b'N', b'H', b'C', 0x01]);
        let actual: Vec<_> = data.iter().collect::<io::Result<_>>()?;

        assert_eq!(actual.len(), 1);

        let (actual_tag, actual_value) = &actual[0];
        assert_eq!(actual_tag, &Tag::ALIGNMENT_HIT_COUNT);
        assert!(matches!(actual_value, Value::UInt8(1)));

        Ok(())
    }
}
