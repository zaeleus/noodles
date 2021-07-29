use std::{convert::TryFrom, ops::Deref, str::FromStr};

use indexmap::IndexSet;
use noodles_vcf::{
    self as vcf,
    header::{Filter, Format, Info, ParseError, Record},
};

/// An indexed map of VCF strings.
///
/// This is also called a dictionary of strings.
///
/// See § 6.2.1 Dictionary of strings (2021-05-13).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StringMap(IndexSet<String>);

impl StringMap {
    fn insert(&mut self, value: String) {
        self.0.insert(value);
    }
}

impl Default for StringMap {
    fn default() -> Self {
        // § 6.2.1 Dictionary of strings (2021-01-13): "Note that 'PASS' is always implicitly
        // encoded as the first entry in the header dictionary."
        Self(vec![Filter::pass().id().into()].into_iter().collect())
    }
}

impl Deref for StringMap {
    type Target = IndexSet<String>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl FromStr for StringMap {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use vcf::header::record::Key;

        let mut string_map = StringMap::default();

        for line in s.lines() {
            if line.starts_with("#CHROM") {
                break;
            }

            let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;

            match record.key() {
                Key::Filter => {
                    let filter = Filter::try_from(record).map_err(ParseError::InvalidFilter)?;
                    string_map.insert(filter.id().into());
                }
                Key::Format => {
                    let format = Format::try_from(record).map_err(ParseError::InvalidFormat)?;
                    string_map.insert(format.id().as_ref().into());
                }
                Key::Info => {
                    let info = Info::try_from(record).map_err(ParseError::InvalidInfo)?;
                    string_map.insert(info.id().as_ref().into());
                }
                _ => {}
            }
        }

        Ok(string_map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        assert_eq!(
            StringMap::default(),
            StringMap(vec![String::from("PASS")].into_iter().collect())
        );
    }

    #[test]
    fn test_from_str() {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20210412
##contig=<ID=sq0,length=8>
##contig=<ID=sq1,length=13>
##contig=<ID=sq2,length=21>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data",IDX=1>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=2>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Quality below 10",IDX=3>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype",IDX=4>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=2>
##ALT=<ID=DEL,Description="Deletion">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse(),
            Ok(StringMap(
                vec![
                    String::from("PASS"),
                    String::from("DP"),
                    String::from("NS"),
                    String::from("q10"),
                    String::from("GT"),
                ]
                .into_iter()
                .collect()
            ))
        );
    }
}
