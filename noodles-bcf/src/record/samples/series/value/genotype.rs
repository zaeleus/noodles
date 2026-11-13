use std::{fmt::Debug, io, iter};

use noodles_vcf::{
    self as vcf, header::FileFormat, variant::record::samples::series::value::genotype::Phasing,
};

use crate::record::{codec::value::Int8, value::array};

const VCF_4_4: FileFormat = FileFormat::new(4, 4);

struct Values<'a>(array::Values<'a, i8>);

impl Values<'_> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<i32>> + '_> {
        Box::new(self.0.iter().map_while(|value| match value {
            Int8::Value(n) => Some(Ok(i32::from(n))),
            Int8::EndOfVector => None,
            Int8::Missing | Int8::Reserved(_) => {
                Some(Err(io::Error::from(io::ErrorKind::InvalidData)))
            }
        }))
    }
}

/// A BCF record samples series genotype value.
pub struct Genotype<'a> {
    file_format: FileFormat,
    values: Values<'a>,
}

impl<'a> Genotype<'a> {
    pub(crate) fn new(file_format: FileFormat, src: &'a [u8]) -> Self {
        Self {
            file_format,
            values: Values(array::Values::new(src)),
        }
    }

    /// Returns an iterator over allele position-phasing pairs.
    pub fn iter(&self) -> impl Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_ {
        let first_allele_phasing = match first_allele_phasing(self.file_format, &self.values) {
            Ok(phasing) => phasing,
            Err(e) => {
                return Box::new(iter::once(Err(e))) as Box<dyn Iterator<Item = io::Result<_>>>;
            }
        };

        Box::new(self.values.iter().enumerate().map(move |(i, result)| {
            let n = result?;
            let position = allele_position(n)?;

            let phasing = if i == 0 {
                first_allele_phasing
            } else {
                allele_phasing(n)
            };

            Ok((position, phasing))
        }))
    }
}

impl Debug for Genotype<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl vcf::variant::record::samples::series::value::Genotype for Genotype<'_> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_> {
        Box::new(self.iter())
    }
}

fn first_allele_phasing(file_format: FileFormat, values: &Values<'_>) -> io::Result<Phasing> {
    // § 6.3.3.9 "Type encoding: Genotype (GT) field" (2024-06-28): "When processing VCF version
    // 4.3 or earlier files, the phasing of the first allele should be treated as missing and
    // inferred from the remaining values."
    if file_format < VCF_4_4 {
        implicit_first_allele_phasing(values)
    } else {
        explicit_first_allele_phasing(values)
    }
}

fn implicit_first_allele_phasing(values: &Values<'_>) -> io::Result<Phasing> {
    for result in values.iter().skip(1) {
        let n = result?;

        if !is_phased(n) {
            return Ok(Phasing::Unphased);
        }
    }

    Ok(Phasing::Phased)
}

fn explicit_first_allele_phasing(values: &Values<'_>) -> io::Result<Phasing> {
    match values.iter().next() {
        Some(result) => result.map(allele_phasing),
        None => Err(io::Error::from(io::ErrorKind::InvalidData)),
    }
}

fn allele_phasing(n: i32) -> Phasing {
    if is_phased(n) {
        Phasing::Phased
    } else {
        Phasing::Unphased
    }
}

fn is_phased(n: i32) -> bool {
    n & 0x01 == 1
}

fn allele_position(n: i32) -> io::Result<Option<usize>> {
    const MISSING: i32 = 0;

    if n < 0 {
        return Err(io::Error::from(io::ErrorKind::InvalidData));
    }

    match n >> 1 {
        MISSING => Ok(None),
        m => usize::try_from(m - 1)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() -> io::Result<()> {
        fn t(src: &[u8], expected: &[(Option<usize>, Phasing)]) -> io::Result<()> {
            let genotype = Genotype::new(FileFormat::default(), src);
            assert_eq!(genotype.iter().collect::<io::Result<Vec<_>>>()?, expected);
            Ok(())
        }

        t(
            &[0x02, 0x02],
            &[(Some(0), Phasing::Unphased), (Some(0), Phasing::Unphased)],
        )?;
        t(
            &[0x02, 0x04],
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Unphased)],
        )?;
        t(
            &[0x04, 0x04],
            &[(Some(1), Phasing::Unphased), (Some(1), Phasing::Unphased)],
        )?;
        t(
            &[0x02, 0x05],
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Phased)],
        )?;
        t(
            &[0x00, 0x00],
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
        )?;
        t(&[0x02], &[(Some(0), Phasing::Unphased)])?;
        t(&[0x04], &[(Some(1), Phasing::Unphased)])?;
        t(
            &[0x02, 0x04, 0x06],
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Unphased),
            ],
        )?;
        t(
            &[0x02, 0x04, 0x07],
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Phased),
            ],
        )?;
        t(&[0x02, 0x81], &[(Some(0), Phasing::Unphased)])?;

        Ok(())
    }

    #[test]
    fn test_first_allele_phasing() -> io::Result<()> {
        const VCF_4_3: FileFormat = FileFormat::new(4, 3);

        let src = [0x02, 0x05]; // |0/1
        let values = Values(array::Values::new(&src));
        assert_eq!(first_allele_phasing(VCF_4_3, &values)?, Phasing::Phased);
        assert_eq!(first_allele_phasing(VCF_4_4, &values)?, Phasing::Unphased);

        let src = [0x03, 0x04]; // /0|1
        let values = Values(array::Values::new(&src));
        assert_eq!(first_allele_phasing(VCF_4_3, &values)?, Phasing::Unphased);
        assert_eq!(first_allele_phasing(VCF_4_4, &values)?, Phasing::Phased);

        Ok(())
    }

    #[test]
    fn test_is_phased() {
        assert!(!is_phased(0x00));
        assert!(is_phased(0x01));
    }

    #[test]
    fn test_allele_position() -> io::Result<()> {
        assert!(allele_position(0x00)?.is_none());
        assert!(allele_position(0x01)?.is_none());
        assert_eq!(allele_position(0x02)?, Some(0));
        assert_eq!(allele_position(0x03)?, Some(0));

        assert!(matches!(allele_position(-1), Err(e) if e.kind() == io::ErrorKind::InvalidData));

        Ok(())
    }
}
