use std::io;

use noodles_vcf as vcf;

/// BCF record alternate bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlternateBases<'a> {
    src: &'a [u8],
    len: usize,
}

impl<'a> AlternateBases<'a> {
    pub(super) fn new(src: &'a [u8], len: usize) -> Self {
        Self { src, len }
    }
}

impl AsRef<[u8]> for AlternateBases<'_> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}

impl vcf::variant::record::AlternateBases for AlternateBases<'_> {
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn len(&self) -> usize {
        self.len
    }

    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<Option<&str>>> + '_> {
        use super::{Value, value::read_value};

        const MISSING: &str = ".";

        let mut src = self.src;

        Box::new((0..self.len()).map(move |_| match read_value(&mut src)? {
            Some(Value::String(None | Some(MISSING))) => Ok(None),
            Some(Value::String(Some(value))) => Ok(Some(value)),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid alt value",
            )),
        }))
    }
}

#[cfg(test)]
mod tests {
    use noodles_vcf::variant::record::AlternateBases as _;

    use super::*;

    #[test]
    fn test_iter() -> io::Result<()> {
        fn t(src: &[u8], len: usize, expected: &[Option<&str>]) -> io::Result<()> {
            let alternate_bases = AlternateBases::new(src, len);
            let actual: Vec<_> = alternate_bases.iter().collect::<io::Result<_>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[], 0, &[])?;
        t(&[0x17, b'C'], 1, &[Some("C")])?;
        t(&[0x17, b'C', 0x07], 2, &[Some("C"), None])?;
        t(&[0x17, b'C', 0x17, b'.'], 2, &[Some("C"), None])?;

        let alternate_bases = AlternateBases::new(&[0x11, 0x00], 1);
        assert!(matches!(
            alternate_bases.iter().collect::<io::Result<Vec<_>>>(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
