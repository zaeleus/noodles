use std::{io, iter};

use crate::variant::record::samples::series::value::genotype::Phasing;

/// VCF record samples series genotype value.
#[derive(Debug)]
pub struct Genotype<'a>(&'a str);

impl<'a> Genotype<'a> {
    pub(crate) fn new(src: &'a str) -> Self {
        Self(src)
    }
}

impl crate::variant::record::samples::series::value::Genotype for Genotype<'_> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_> {
        let mut src = self.0;

        let first_allele_result = parse_first_allele(&mut src);

        Box::new(
            iter::once(first_allele_result).chain(iter::from_fn(move || {
                if src.is_empty() {
                    None
                } else {
                    Some(parse_allele(&mut src))
                }
            })),
        )
    }
}

fn parse_first_allele(src: &mut &str) -> io::Result<(Option<usize>, Phasing)> {
    let mut buf = next_allele(src);

    let phasing = if buf.starts_with(['|', '/']) {
        let p = parse_phasing(&buf[..1])?;
        buf = &buf[1..];
        p
    } else if src
        .as_bytes()
        .iter()
        .filter(|&&b| matches!(b, b'|' | b'/'))
        .any(|&b| b == b'/')
    {
        Phasing::Unphased
    } else {
        Phasing::Phased
    };

    let position = parse_position(buf)?;

    Ok((position, phasing))
}

fn parse_allele(src: &mut &str) -> io::Result<(Option<usize>, Phasing)> {
    let buf = next_allele(src);

    let phasing = parse_phasing(&buf[..1])?;
    let position = parse_position(&buf[1..])?;

    Ok((position, phasing))
}

fn next_allele<'a>(src: &mut &'a str) -> &'a str {
    let (buf, rest) = match src.chars().skip(1).position(is_phasing_indicator) {
        Some(i) => src.split_at(i + 1),
        None => src.split_at(src.len()),
    };

    *src = rest;

    buf
}

fn is_phasing_indicator(c: char) -> bool {
    const PHASED: char = '|';
    const UNPHASED: char = '/';

    matches!(c, PHASED | UNPHASED)
}

fn parse_phasing(src: &str) -> io::Result<Phasing> {
    const PHASED: &str = "|";
    const UNPHASED: &str = "/";

    match src {
        PHASED => Ok(Phasing::Phased),
        UNPHASED => Ok(Phasing::Unphased),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid phasing indicator",
        )),
    }
}

fn parse_position(src: &str) -> io::Result<Option<usize>> {
    const MISSING: &str = ".";

    match src {
        MISSING => Ok(None),
        _ => src
            .parse()
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::variant::record::samples::series::value::Genotype as _;

    #[test]
    fn test_iter() -> io::Result<()> {
        fn t(src: &str, expected: &[(Option<usize>, Phasing)]) -> io::Result<()> {
            let genotype = Genotype::new(src);
            let actual: Vec<_> = genotype.iter().collect::<io::Result<_>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(
            "0|0",
            &[(Some(0), Phasing::Phased), (Some(0), Phasing::Phased)],
        )?;
        t(
            "0/1",
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Unphased)],
        )?;
        t(
            "|0/1",
            &[(Some(0), Phasing::Phased), (Some(1), Phasing::Unphased)],
        )?;
        t(
            "|1/2|3",
            &[
                (Some(1), Phasing::Phased),
                (Some(2), Phasing::Unphased),
                (Some(3), Phasing::Phased),
            ],
        )?;
        t(
            "./.",
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
        )?;

        Ok(())
    }
}
