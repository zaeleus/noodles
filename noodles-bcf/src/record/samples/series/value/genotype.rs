use std::{fmt::Debug, io};

use noodles_vcf::{
    self as vcf, header::FileFormat, variant::record::samples::series::value::genotype::Phasing,
};

use crate::record::codec::value::Int8;

const VCF_4_4: FileFormat = FileFormat::new(4, 4);

/// A BCF record samples series genotype value.
pub struct Genotype<'a> {
    file_format: FileFormat,
    src: &'a [u8],
}

impl<'a> Genotype<'a> {
    pub(crate) fn new(file_format: FileFormat, src: &'a [u8]) -> Self {
        Self { file_format, src }
    }

    /// Returns an iterator over allele position-phasing pairs.
    pub fn iter(&self) -> impl Iterator<Item = (Option<usize>, Phasing)> + '_ {
        const MISSING: i8 = -1;

        let first_allele_phasing = first_allele_phasing(self.file_format, self.src);

        self.src
            .iter()
            .map(|n| *n as i8)
            .enumerate()
            .take_while(|(_, n)| matches!(Int8::from(*n), Int8::Value(_)))
            .map(move |(i, n)| {
                let j = (n >> 1) - 1;

                let position = if j == MISSING { None } else { Some(j as usize) };

                let phasing = if i == 0 {
                    first_allele_phasing
                } else {
                    allele_phasing(n)
                };

                (position, phasing)
            })
    }
}

impl Debug for Genotype<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl vcf::variant::record::samples::series::value::Genotype for Genotype<'_> {
    fn iter(&self) -> Box<dyn Iterator<Item = io::Result<(Option<usize>, Phasing)>> + '_> {
        Box::new(self.iter().map(Ok))
    }
}

fn first_allele_phasing(file_format: FileFormat, src: &[u8]) -> Phasing {
    // § 6.3.3.9 "Type encoding: Genotype (GT) field" (2024-06-28): "When processing VCF version
    // 4.3 or earlier files, the phasing of the first allele should be treated as missing and
    // inferred from the remaining values."
    if file_format < VCF_4_4 {
        implicit_first_allele_phasing(src)
    } else {
        explicit_first_allele_phasing(src)
    }
}

fn implicit_first_allele_phasing(src: &[u8]) -> Phasing {
    let mut phasing = Phasing::Phased;

    for &n in src.iter().skip(1) {
        let n = n as i8;

        if !matches!(Int8::from(n), Int8::Value(_)) {
            break;
        }

        if !is_phased(n) {
            phasing = Phasing::Unphased;
            break;
        }
    }

    phasing
}

fn explicit_first_allele_phasing(src: &[u8]) -> Phasing {
    allele_phasing(src[0] as i8)
}

fn allele_phasing(n: i8) -> Phasing {
    if is_phased(n) {
        Phasing::Phased
    } else {
        Phasing::Unphased
    }
}

fn is_phased(n: i8) -> bool {
    n & 0x01 == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        fn t(src: &[u8], expected: &[(Option<usize>, Phasing)]) {
            let genotype = Genotype::new(FileFormat::default(), src);
            assert_eq!(genotype.iter().collect::<Vec<_>>(), expected);
        }

        t(
            &[0x02, 0x02],
            &[(Some(0), Phasing::Unphased), (Some(0), Phasing::Unphased)],
        );
        t(
            &[0x02, 0x04],
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Unphased)],
        );
        t(
            &[0x04, 0x04],
            &[(Some(1), Phasing::Unphased), (Some(1), Phasing::Unphased)],
        );
        t(
            &[0x02, 0x05],
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Phased)],
        );
        t(
            &[0x00, 0x00],
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
        );
        t(&[0x02], &[(Some(0), Phasing::Unphased)]);
        t(&[0x04], &[(Some(1), Phasing::Unphased)]);
        t(
            &[0x02, 0x04, 0x06],
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Unphased),
            ],
        );
        t(
            &[0x02, 0x04, 0x07],
            &[
                (Some(0), Phasing::Unphased),
                (Some(1), Phasing::Unphased),
                (Some(2), Phasing::Phased),
            ],
        );
        t(&[0x02, 0x81], &[(Some(0), Phasing::Unphased)]);
    }

    #[test]
    fn test_first_allele_phasing() {
        const VCF_4_3: FileFormat = FileFormat::new(4, 3);

        let src = [0x02, 0x05]; // |0/1
        assert_eq!(first_allele_phasing(VCF_4_3, &src), Phasing::Phased);
        assert_eq!(first_allele_phasing(VCF_4_4, &src), Phasing::Unphased);

        let src = [0x03, 0x04]; // /0|1
        assert_eq!(first_allele_phasing(VCF_4_3, &src), Phasing::Unphased);
        assert_eq!(first_allele_phasing(VCF_4_4, &src), Phasing::Phased);
    }

    #[test]
    fn test_is_phased() {
        assert!(!is_phased(0x00));
        assert!(is_phased(0x01));
    }
}
