use std::{fmt::Debug, io};

use noodles_vcf::{self as vcf, variant::record::samples::series::value::genotype::Phasing};

use crate::record::codec::value::Int8;

/// A BCF record samples series genotype value.
pub struct Genotype<'a>(&'a [u8]);

impl<'a> Genotype<'a> {
    pub(crate) fn new(src: &'a [u8]) -> Self {
        Self(src)
    }

    /// Returns an iterator over allele position-phasing pairs.
    pub fn iter(&self) -> impl Iterator<Item = (Option<usize>, Phasing)> + '_ {
        const MISSING: i8 = -1;

        let mut first_allele_phasing = Phasing::Phased;

        for &n in self.0.iter().skip(1) {
            let n = n as i8;

            if !matches!(Int8::from(n), Int8::Value(_)) {
                break;
            }

            if !is_phased(n) {
                first_allele_phasing = Phasing::Unphased;
                break;
            }
        }

        self.0
            .iter()
            .map(|n| *n as i8)
            .enumerate()
            .take_while(|(_, n)| matches!(Int8::from(*n), Int8::Value(_)))
            .map(move |(i, n)| {
                let j = (n >> 1) - 1;

                let position = if j == MISSING { None } else { Some(j as usize) };

                let phasing = if i == 0 {
                    first_allele_phasing
                } else if is_phased(n) {
                    Phasing::Phased
                } else {
                    Phasing::Unphased
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

fn is_phased(n: i8) -> bool {
    n & 0x01 == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iter() {
        fn t(src: &[u8], expected: &[(Option<usize>, Phasing)]) {
            let genotype = Genotype::new(src);
            assert_eq!(genotype.iter().collect::<Vec<_>>(), expected);
        }

        // § 6.3.3 "Type encoding" (2023-08-23)
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
            &[(Some(0), Phasing::Phased), (Some(1), Phasing::Phased)],
        );
        t(
            &[0x00, 0x00],
            &[(None, Phasing::Unphased), (None, Phasing::Unphased)],
        );
        t(&[0x02], &[(Some(0), Phasing::Phased)]);
        t(&[0x04], &[(Some(1), Phasing::Phased)]);
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
        t(&[0x02, 0x81], &[(Some(0), Phasing::Phased)]);
    }

    #[test]
    fn test_is_phased() {
        assert!(!is_phased(0x00));
        assert!(is_phased(0x01));
    }
}
