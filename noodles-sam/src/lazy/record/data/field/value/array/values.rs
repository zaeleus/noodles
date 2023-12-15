use std::{io, marker::PhantomData};

use lexical_core::FromLexical;

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

    /// Returns an iterator over values.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<N>> + '_ {
        self.src.split(delimiter).map(parse_num)
    }
}

impl<'a, N> crate::alignment::record::data::field::value::array::Values<'a, N> for Values<'a, N>
where
    N: FromLexical,
{
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<N>> + '_> {
        Box::new(self.iter())
    }
}

fn delimiter(b: &u8) -> bool {
    const DELIMITER: u8 = b',';
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
    fn test_iter() -> io::Result<()> {
        let values = Values::<'_, u8>::new(b"8,13");
        let actual: Vec<_> = values.iter().collect::<Result<_, _>>()?;
        assert_eq!(actual, [8, 13]);
        Ok(())
    }
}
