//! An indexed map of VCF strings.

mod string_map;

use std::str::{FromStr, Lines};

pub use self::string_map::StringMap;
use crate::{
    header::{
        parser::{parse_record, Entry},
        FileFormat, ParseError, Record,
    },
    Header,
};

/// An indexed map of VCF strings (FILTER, FORMAT, and INFO).
pub type StringStringMap = StringMap;

/// An indexed map of VCF contig names.
pub type ContigStringMap = StringMap;

/// An indexed map of VCF strings.
///
/// This includes both the dictionary of strings and dictionary of contigs.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct StringMaps {
    string_string_map: StringStringMap,
    contig_string_map: ContigStringMap,
}

impl StringMaps {
    /// Returns an indexed map of VCF strings (FILTER, FORMAT, and INFO).
    ///
    /// The filter ID "PASS" is always the first entry in the string string map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{
    ///         record::value::{map::{Contig, Filter, Format, Info}, Map},
    ///         StringMaps,
    ///     },
    ///     variant::{record::info, record_buf::samples},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(info::field::key::TOTAL_DEPTH, Map::<Info>::from(info::field::key::TOTAL_DEPTH))
    ///     .add_filter("q10", Map::<Filter>::new("Quality below 10"))
    ///     .add_format(samples::keys::key::READ_DEPTH, Map::<Format>::from(samples::keys::key::READ_DEPTH))
    ///     .add_contig("sq0", Map::<Contig>::new())
    ///     .build();
    ///
    /// let string_maps = StringMaps::try_from(&header)?;
    /// let string_string_map = string_maps.strings();
    ///
    /// assert_eq!(string_string_map.get_index(0), Some("PASS"));
    /// assert_eq!(string_string_map.get_index(1), Some("DP"));
    /// assert_eq!(string_string_map.get_index(2), Some("q10"));
    /// assert!(string_string_map.get_index(3).is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn strings(&self) -> &StringStringMap {
        &self.string_string_map
    }

    fn strings_mut(&mut self) -> &mut StringStringMap {
        &mut self.string_string_map
    }

    /// Returns an indexed map of contig names.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::{
    ///         record::value::{map::{Contig, Filter, Format, Info}, Map},
    ///         StringMaps,
    ///     },
    ///     variant::{record::info, record_buf::samples},
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(info::field::key::TOTAL_DEPTH, Map::<Info>::from(info::field::key::TOTAL_DEPTH))
    ///     .add_filter("q10", Map::<Filter>::new("Quality below 10"))
    ///     .add_format(samples::keys::key::READ_DEPTH, Map::<Format>::from(samples::keys::key::READ_DEPTH))
    ///     .add_contig("sq0", Map::<Contig>::new())
    ///     .build();
    ///
    /// let string_maps = StringMaps::try_from(&header)?;
    /// let contig_string_map = string_maps.contigs();
    ///
    /// assert_eq!(contig_string_map.get_index(0), Some("sq0"));
    /// assert!(contig_string_map.get_index(1).is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn contigs(&self) -> &ContigStringMap {
        &self.contig_string_map
    }

    fn contigs_mut(&mut self) -> &mut ContigStringMap {
        &mut self.contig_string_map
    }

    #[doc(hidden)]
    pub fn insert_entry(&mut self, entry: &Entry<'_>) -> Result<(), ParseError> {
        match entry {
            Entry::Contig(id, contig) => insert(self.contigs_mut(), id, contig.idx()),
            Entry::Filter(id, filter) => insert(self.strings_mut(), id, filter.idx()),
            Entry::Format(id, format) => insert(self.strings_mut(), id, format.idx()),
            Entry::Info(id, info) => insert(self.strings_mut(), id, info.idx()),
            _ => Ok(()),
        }
    }
}

impl Default for StringMaps {
    fn default() -> Self {
        // § 6.2.1 Dictionary of strings (2021-01-13): "Note that 'PASS' is always implicitly
        // encoded as the first entry in the header dictionary."
        let mut string_string_map = StringMap::default();
        string_string_map.insert(String::from("PASS"));

        let contig_string_map = StringMap::default();

        Self {
            string_string_map,
            contig_string_map,
        }
    }
}

impl FromStr for StringMaps {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut string_maps = Self::default();

        let mut lines = s.lines();
        let file_format = parse_file_format(&mut lines)?;

        for line in &mut lines {
            if line.starts_with("#CHROM") {
                break;
            }

            let record =
                parse_record(line.as_bytes(), file_format).map_err(ParseError::InvalidRecord)?;

            match record {
                Record::Contig(id, contig) => {
                    insert(string_maps.contigs_mut(), id.as_ref(), contig.idx())?;
                }
                Record::Filter(id, filter) => {
                    insert(string_maps.strings_mut(), &id, filter.idx())?;
                }
                Record::Format(id, format) => {
                    insert(string_maps.strings_mut(), id.as_ref(), format.idx())?;
                }
                Record::Info(id, info) => {
                    insert(string_maps.strings_mut(), id.as_ref(), info.idx())?;
                }
                _ => {}
            }
        }

        Ok(string_maps)
    }
}

fn parse_file_format(lines: &mut Lines<'_>) -> Result<FileFormat, ParseError> {
    let record = lines
        .next()
        .ok_or(ParseError::MissingFileFormat)
        .and_then(|line| {
            parse_record(line.as_bytes(), Default::default()).map_err(ParseError::InvalidRecord)
        })?;

    match record {
        Record::FileFormat(file_format) => Ok(file_format),
        _ => Err(ParseError::MissingFileFormat),
    }
}

fn insert(string_map: &mut StringMap, id: &str, idx: Option<usize>) -> Result<(), ParseError> {
    if let Some(i) = idx {
        if let Some((j, entry)) = string_map.get_full(id) {
            let actual = (i, id.into());
            let expected = (j, entry.into());

            if actual != expected {
                return Err(ParseError::StringMapPositionMismatch(actual, expected));
            }
        } else {
            string_map.insert_at(i, id.into());
        }
    } else {
        string_map.insert(id.into());
    }

    Ok(())
}

impl TryFrom<&Header> for StringMaps {
    type Error = ParseError;

    fn try_from(header: &Header) -> Result<Self, Self::Error> {
        let mut string_maps = StringMaps::default();

        for (id, contig) in header.contigs() {
            insert(string_maps.contigs_mut(), id.as_ref(), contig.idx())?;
        }

        for (id, info) in header.infos() {
            insert(string_maps.strings_mut(), id.as_ref(), info.idx())?;
        }

        for (id, filter) in header.filters() {
            insert(string_maps.strings_mut(), id, filter.idx())?;
        }

        for (id, format) in header.formats() {
            insert(string_maps.strings_mut(), id.as_ref(), format.idx())?;
        }

        Ok(string_maps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let actual = StringMaps::default();

        let mut string_string_map = StringMap::default();
        string_string_map.insert("PASS".into());

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(actual, expected);
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

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("DP"), 2),
                (String::from("q10"), 3),
                (String::from("GT"), 4),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                Some(String::from("DP")),
                Some(String::from("q10")),
                Some(String::from("GT")),
            ],
        };

        let contig_string_map = StringMap {
            indices: [
                (String::from("sq0"), 0),
                (String::from("sq1"), 1),
                (String::from("sq2"), 2),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("sq0")),
                Some(String::from("sq1")),
                Some(String::from("sq2")),
            ],
        };

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_file_format() {
        // FORMAT MQ is an `Integer` in VCF 4.2 and `Float` in VCF 4.3.
        let s = r#"##fileformat=VCFv4.2
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="RMS mapping quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let mut string_string_map = StringMap::default();
        string_string_map.insert(String::from("PASS"));
        string_string_map.insert(String::from("MQ"));

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_mixed_positions() {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20210412
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data",IDX=1>
##FILTER=<ID=PASS,Description="All filters passed",IDX=0>
##FILTER=<ID=q10,Description="Quality below 10",IDX=3>
##FILTER=<ID=q15,Description="Quality below 15",IDX=4>
##FILTER=<ID=q20,Description="Quality below 20">
##FILTER=<ID=NS,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("q10"), 3),
                (String::from("q15"), 4),
                (String::from("q20"), 5),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                None,
                Some(String::from("q10")),
                Some(String::from("q15")),
                Some(String::from("q20")),
            ],
        };

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(s.parse(), Ok(expected));
    }

    #[test]
    fn test_from_str_with_a_position_mismatch() {
        let s = r#"##fileformat=VCFv4.3
##FILTER=<ID=PASS,Description="All filters passed",IDX=8>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMaps>(),
            Err(ParseError::StringMapPositionMismatch(
                (8, String::from("PASS")),
                (0, String::from("PASS"))
            ))
        );

        let s = r#"##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Combined depth across samples",IDX=1>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth",IDX=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        assert_eq!(
            s.parse::<StringMaps>(),
            Err(ParseError::StringMapPositionMismatch(
                (2, String::from("DP")),
                (1, String::from("DP"))
            ))
        );
    }

    #[test]
    fn test_try_from_vcf_header_for_string_maps() -> Result<(), Box<dyn std::error::Error>> {
        use crate::{
            header::record::value::{
                map::{AlternativeAllele, Contig, Filter, Format, Info},
                Map,
            },
            variant::record::info,
            variant::record_buf::samples,
        };

        let header = Header::builder()
            .add_contig("sq0", Map::<Contig>::new())
            .add_contig("sq1", Map::<Contig>::new())
            .add_contig("sq2", Map::<Contig>::new())
            .add_info(
                info::field::key::SAMPLES_WITH_DATA_COUNT,
                Map::<Info>::from(info::field::key::SAMPLES_WITH_DATA_COUNT),
            )
            .add_info(
                info::field::key::TOTAL_DEPTH,
                Map::<Info>::from(info::field::key::TOTAL_DEPTH),
            )
            .add_filter("PASS", Map::<Filter>::pass())
            .add_filter("q10", Map::<Filter>::new("Quality below 10"))
            .add_format(
                samples::keys::key::GENOTYPE,
                Map::<Format>::from(samples::keys::key::GENOTYPE),
            )
            .add_format(
                samples::keys::key::READ_DEPTH,
                Map::<Format>::from(samples::keys::key::READ_DEPTH),
            )
            .add_alternative_allele("DEL", Map::<AlternativeAllele>::new("Deletion"))
            .build();

        let actual = StringMaps::try_from(&header)?;

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("DP"), 2),
                (String::from("q10"), 3),
                (String::from("GT"), 4),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                Some(String::from("DP")),
                Some(String::from("q10")),
                Some(String::from("GT")),
            ],
        };

        let contig_string_map = StringMap {
            indices: [
                (String::from("sq0"), 0),
                (String::from("sq1"), 1),
                (String::from("sq2"), 2),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("sq0")),
                Some(String::from("sq1")),
                Some(String::from("sq2")),
            ],
        };

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_try_from_vcf_header_for_string_maps_with_idx() -> Result<(), Box<dyn std::error::Error>>
    {
        use crate::{
            header::record::value::{
                map::{Filter, Info},
                Map,
            },
            variant::record::info,
        };

        let ns = {
            let mut map = Map::<Info>::from(info::field::key::SAMPLES_WITH_DATA_COUNT);
            *map.idx_mut() = Some(1);
            map
        };

        let header = Header::builder()
            .add_filter(
                "PASS",
                Map::<Filter>::builder()
                    .set_description("All filters passed")
                    .set_idx(0)
                    .build()?,
            )
            .add_filter(
                "q10",
                Map::<Filter>::builder()
                    .set_description("Quality below 10")
                    .set_idx(3)
                    .build()?,
            )
            .add_filter(
                "q15",
                Map::<Filter>::builder()
                    .set_description("Quality below 15")
                    .set_idx(4)
                    .build()?,
            )
            .add_filter("q20", Map::<Filter>::new("Quality below 20"))
            .add_filter("NS", Map::<Filter>::new(""))
            .add_info(info::field::key::SAMPLES_WITH_DATA_COUNT, ns)
            .build();

        let actual = StringMaps::try_from(&header)?;

        let string_string_map = StringMap {
            indices: [
                (String::from("PASS"), 0),
                (String::from("NS"), 1),
                (String::from("q10"), 3),
                (String::from("q15"), 4),
                (String::from("q20"), 5),
            ]
            .into_iter()
            .collect(),
            entries: vec![
                Some(String::from("PASS")),
                Some(String::from("NS")),
                None,
                Some(String::from("q10")),
                Some(String::from("q15")),
                Some(String::from("q20")),
            ],
        };

        let contig_string_map = StringMap::default();

        let expected = StringMaps {
            string_string_map,
            contig_string_map,
        };

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_file_format() {
        let s = "##fileformat=VCFv4.3\n";
        let mut lines = s.lines();
        assert_eq!(parse_file_format(&mut lines), Ok(FileFormat::new(4, 3)));

        let s = "";
        let mut lines = s.lines();
        assert_eq!(
            parse_file_format(&mut lines),
            Err(ParseError::MissingFileFormat)
        );

        let s = "##fileDate=20211119";
        let mut lines = s.lines();
        assert_eq!(
            parse_file_format(&mut lines),
            Err(ParseError::MissingFileFormat)
        );

        let s = "fileformat=VCFv4.3";
        let mut lines = s.lines();
        assert!(matches!(
            parse_file_format(&mut lines),
            Err(ParseError::InvalidRecord(_))
        ));

        let s = "##fileformat=VCF43\n";
        let mut lines = s.lines();
        assert!(matches!(
            parse_file_format(&mut lines),
            Err(ParseError::InvalidRecord(_))
        ));
    }
}
