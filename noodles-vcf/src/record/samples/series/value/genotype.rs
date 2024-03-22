use std::{io, iter};

/// VCF record samples series genotype value.
#[derive(Debug)]
pub struct Genotype<'a>(&'a str);

impl<'a> Genotype<'a> {
    pub(crate) fn new(src: &'a str) -> Self {
        Self(src)
    }
}

impl<'a> crate::variant::record::samples::series::value::Genotype for Genotype<'a> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, u8)>> + '_> {
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

fn parse_first_allele(src: &mut &str) -> io::Result<(Option<usize>, u8)> {
    let mut buf = next_allele(src);

    let phasing = if buf.starts_with(|c| matches!(c, '|' | '/')) {
        let p = buf.as_bytes()[0];
        buf = &buf[1..];
        p
    } else if src
        .as_bytes()
        .iter()
        .filter(|&&b| matches!(b, b'|' | b'/'))
        .any(|&b| b == b'/')
    {
        b'/'
    } else {
        b'|'
    };

    let position = parse_position(buf)?;

    Ok((position, phasing))
}

fn parse_allele(src: &mut &str) -> io::Result<(Option<usize>, u8)> {
    let buf = next_allele(src);

    let phasing = buf.as_bytes()[0];
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
        fn t(src: &str, expected: &[(Option<usize>, u8)]) -> io::Result<()> {
            let genotype = Genotype::new(src);
            let actual: Vec<_> = genotype.iter().collect::<io::Result<_>>()?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t("0|0", &[(Some(0), b'|'), (Some(0), b'|')])?;
        t("0/1", &[(Some(0), b'/'), (Some(1), b'/')])?;
        t("|0/1", &[(Some(0), b'|'), (Some(1), b'/')])?;
        t(
            "|1/2|3",
            &[(Some(1), b'|'), (Some(2), b'/'), (Some(3), b'|')],
        )?;
        t("./.", &[(None, b'/'), (None, b'/')])?;

        Ok(())
    }
}
