use std::{io, iter, marker::PhantomData};

use lexical_core::FromLexical;

const DELIMITER: u8 = b',';

#[derive(Debug, PartialEq)]
pub struct Values<'a, N> {
    src: &'a [u8],
    _marker: PhantomData<N>,
}

impl<'a, N> Values<'a, N>
where
    N: FromLexical,
{
    pub(crate) fn new(src: &'a [u8]) -> Self {
        Self {
            src,
            _marker: PhantomData,
        }
    }

    fn len(&self) -> usize {
        if self.src.is_empty() {
            0
        } else {
            self.src.iter().filter(|&&b| b == DELIMITER).count() + 1
        }
    }

    /// Returns an iterator over values.
    pub fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        if self.src.is_empty() {
            Box::new(iter::empty())
        } else {
            Box::new(self.src.split(delimiter).map(parse_num))
        }
    }
}

impl<'a, N> crate::alignment::record::data::field::value::array::Values<'a, N> for Values<'a, N>
where
    N: FromLexical,
{
    fn len(&self) -> usize {
        self.len()
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        self.iter()
    }
}

fn delimiter(b: &u8) -> bool {
    *b == DELIMITER
}

fn parse_num<N>(src: &[u8]) -> io::Result<N>
where
    N: FromLexical,
{
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_len() {
        assert_eq!(Values::<'_, u8>::new(b"").len(), 0);
        assert_eq!(Values::<'_, u8>::new(b"0").len(), 1);
        assert_eq!(Values::<'_, u8>::new(b"0,0").len(), 2);
    }

    #[test]
    fn test_iter() -> io::Result<()> {
        assert!(Values::<'_, u8>::new(b"").iter().next().is_none());
        assert_eq!(
            Values::<'_, u8>::new(b"5")
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5]
        );
        assert_eq!(
            Values::<'_, u8>::new(b"5,8")
                .iter()
                .collect::<io::Result<Vec<_>>>()?,
            [5, 8]
        );

        Ok(())
    }
}
