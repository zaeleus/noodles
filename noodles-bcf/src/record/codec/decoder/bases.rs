use std::{io, num::NonZero};

use noodles_vcf::variant::record_buf::AlternateBases;

use super::read_value;
use crate::record::codec::Value;

pub(crate) fn read_ref_alt(
    src: &mut &[u8],
    len: NonZero<usize>,
) -> io::Result<(String, AlternateBases)> {
    let mut alleles = Vec::with_capacity(len.get());

    for _ in 0..len.get() {
        match read_value(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))? {
            Some(Value::String(Some(s))) => alleles.push(s),
            Some(Value::String(None)) => alleles.push("."),
            v => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid ref_alt: expected string, got {v:?}"),
                ));
            }
        }
    }

    let (raw_reference_bases, raw_alternate_bases) = alleles
        .split_first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing reference bases"))?;

    let reference_bases = raw_reference_bases.to_string();
    let alternate_bases = raw_alternate_bases.iter().map(|s| s.to_string()).collect();

    Ok((reference_bases, alternate_bases))
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
            0x17, 0x43, // alt = C
        ];
        let actual = read_ref_alt(&mut &src[..], const { NonZero::new(2).unwrap() })?;
        let expected = (String::from("A"), [String::from("C")].into_iter().collect());
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
