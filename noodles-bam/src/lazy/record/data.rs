pub mod field;

use std::{borrow::Borrow, io, iter};

use noodles_sam as sam;

use self::field::{decode_field, Tag, Value};

/// Raw BAM record data.
#[derive(Debug, Eq, PartialEq)]
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

impl<'a> AsRef<[u8]> for Data<'a> {
    fn as_ref(&self) -> &[u8] {
        self.0
    }
}

impl<'a> TryFrom<Data<'a>> for sam::record::Data {
    type Error = io::Error;

    fn try_from(bam_data: Data<'a>) -> Result<Self, Self::Error> {
        use crate::record::codec::decoder::get_data;

        let mut src = bam_data.0;
        let mut sam_data = Self::default();
        get_data(&mut src, &mut sam_data)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        Ok(sam_data)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get() -> io::Result<()> {
        use sam::record::data::field::tag;

        let data = Data::new(&[b'N', b'H', b'C', 0x01]);

        assert!(data.get(&tag::ALIGNMENT_HIT_COUNT).is_some());
        assert!(data.get(&[b'N', b'H']).is_some());

        assert!(data.get(&tag::COMMENT).is_none());

        Ok(())
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        let data = Data::new(&[]);
        assert!(data.iter().next().is_none());

        let data = Data::new(&[b'N', b'H', b'C', 0x01]);
        let actual: Vec<_> = data.iter().collect::<io::Result<_>>()?;
        let expected = [([b'N', b'H'], Value::UInt8(1))];
        assert_eq!(actual, expected);

        Ok(())
    }
}
