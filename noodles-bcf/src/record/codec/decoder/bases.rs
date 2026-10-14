use std::{io, num::NonZero};

use noodles_vcf::variant::record_buf::AlternateBases;

use super::read_value;
use crate::record::codec::Value;

const MISSING: &str = ".";

pub(crate) fn read_ref_alt(
    src: &mut &[u8],
    len: NonZero<usize>,
) -> io::Result<(String, AlternateBases)> {
    let raw_reference_bases = read_allele(src)?.unwrap_or(MISSING);

    let alternate_bases_count = len.get() - 1;

    let alternate_bases = (0..alternate_bases_count)
        .map(|_| read_allele(src).map(|s| s.unwrap_or(MISSING).to_string()))
        .collect::<io::Result<_>>()?;

    Ok((raw_reference_bases.into(), alternate_bases))
}

fn read_allele<'a>(src: &mut &'a [u8]) -> io::Result<Option<&'a str>> {
    let value = read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    match value {
        Some(Value::String(s)) => Ok(s),
        v => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid allele: expected string, got {v:?}"),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_ref_alt() -> io::Result<()> {
        let src = [0x17, 0x41]; // ref = A
        let actual = read_ref_alt(&mut &src[..], NonZero::<usize>::MIN)?;
        let expected = (String::from("A"), AlternateBases::default());
        assert_eq!(actual, expected);

        let src = [
            0x17, 0x41, // ref = A
            0x17, 0x43, 0x07, // alt = C,.
        ];
        let actual = read_ref_alt(&mut &src[..], const { NonZero::new(3).unwrap() })?;
        let expected = (
            String::from("A"),
            [String::from("C"), String::from(".")].into_iter().collect(),
        );
        assert_eq!(actual, expected);

        assert!(matches!(
            read_ref_alt(&mut &[][..], NonZero::<usize>::MIN),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            read_ref_alt(&mut &[0x11, 0x00][..], NonZero::<usize>::MIN),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
