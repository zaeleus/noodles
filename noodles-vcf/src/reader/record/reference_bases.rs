mod base;

use std::io;

use self::base::parse_base;
use crate::record::ReferenceBases;

pub(super) fn parse_reference_bases(
    s: &str,
    reference_bases: &mut ReferenceBases,
) -> io::Result<()> {
    if s.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty reference bases",
        ));
    }

    reference_bases.0.clear();

    for c in s.chars() {
        let base = parse_base(c).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        reference_bases.0.push(base);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_reference_bases() -> Result<(), Box<dyn std::error::Error>> {
        use crate::record::reference_bases::Base;

        let mut reference_bases = ReferenceBases::try_from(vec![Base::N])?;

        let expected = [Base::A, Base::T, Base::C, Base::G, Base::N];

        parse_reference_bases("ATCGN", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        parse_reference_bases("atcgn", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        parse_reference_bases("AtCgN", &mut reference_bases)?;
        assert_eq!(&reference_bases[..], &expected[..]);

        assert!(matches!(
            parse_reference_bases("", &mut reference_bases),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_reference_bases(".", &mut reference_bases),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        assert!(matches!(
            parse_reference_bases("Z", &mut reference_bases),
            Err(e) if e.kind() == io::ErrorKind::InvalidData,
        ));

        Ok(())
    }
}
