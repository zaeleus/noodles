use std::error;

use indexmap::IndexSet;

use super::{
    alternative_allele::{self, AlternativeAllele},
    contig::{self, Contig},
    file_format::{self, FileFormat},
    filter::{self, Filter},
    format::{self, Format},
    info::{self, Info},
    meta::{self, Meta},
    pedigree::{self, Pedigree},
    record::{self, Record},
    sample::{self, Sample},
    Builder, Header,
};

/// An error returned when a raw VCF header fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The file format (`fileformat`) is missing.
    MissingFileFormat,
    /// The file format (`fileformat`) appears other than the first line.
    UnexpectedFileFormat,
    /// The file format (`fileformat`) is invalid.
    InvalidFileFormat(file_format::ParseError),
    /// A record is invalid.
    InvalidRecord(record::ParseError),
    /// A record has an invalid value.
    InvalidRecordValue,
    /// An information record (`INFO`) is invalid.
    InvalidInfo(info::TryFromRecordError),
    /// A filter record (`FILTER`) is invalid.
    InvalidFilter(filter::TryFromRecordError),
    /// A genotype format record (`FORMAT`) is invalid.
    InvalidFormat(format::TryFromRecordError),
    /// A symboloic alternate allele record (`ALT`) is invalid.
    InvalidAlternativeAllele(alternative_allele::TryFromRecordError),
    /// A contig record (`contig`) is invalid.
    InvalidContig(contig::TryFromRecordError),
    /// A meta record (`META`) is invalid.
    InvalidMeta(meta::TryFromRecordError),
    /// A sample record (`SAMPLE`) is invalid.
    InvalidSample(sample::TryFromRecordError),
    /// A pedigree record (`PEDIGREE`) is invalid.
    InvalidPedigree(pedigree::TryFromRecordError),
    /// The header is missing.
    MissingHeader,
    /// The header is invalid.
    InvalidHeader(String, String),
    /// A sample name is duplicated.
    ///
    /// § 1.5 Header line syntax (2021-01-13): "Duplicate sample IDs are not allowed."
    DuplicateSampleName(String),
    /// More data unexpectedly appears after the header header (`#CHROM`...).
    ExpectedEof,
    /// The position of the entry in the string match does not match the absolute position defined
    /// by the `IDX` field of a record.
    StringMapPositionMismatch((usize, String), (usize, String)),
}

impl error::Error for ParseError {}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingFileFormat => f.write_str("missing fileformat"),
            Self::UnexpectedFileFormat => f.write_str("unexpected file format"),
            Self::InvalidFileFormat(e) => write!(f, "invalid file format: {}", e),
            Self::InvalidRecord(e) => write!(f, "invalid record: {}", e),
            Self::InvalidRecordValue => f.write_str("invalid record value"),
            Self::InvalidInfo(e) => write!(f, "invalid info: {}", e),
            Self::InvalidFilter(e) => write!(f, "invalid filter: {}", e),
            Self::InvalidFormat(e) => write!(f, "invalid format: {}", e),
            Self::InvalidAlternativeAllele(e) => {
                write!(f, "invalid alternative allele: {}", e)
            }
            Self::InvalidContig(e) => write!(f, "invalid contig: {}", e),
            Self::InvalidMeta(e) => write!(f, "invalid meta: {}", e),
            Self::InvalidSample(e) => write!(f, "invalid sample: {}", e),
            Self::InvalidPedigree(e) => write!(f, "invalid pedigree: {}", e),
            Self::MissingHeader => f.write_str("missing header"),
            Self::InvalidHeader(actual, expected) => {
                write!(f, "invalid header: expected {}, got {}", expected, actual)
            }
            Self::DuplicateSampleName(sample_name) => {
                write!(f, "duplicate sample name: {}", sample_name)
            }
            Self::ExpectedEof => f.write_str("expected EOF"),
            Self::StringMapPositionMismatch(actual, expected) => write!(
                f,
                "string map position mismatch: expected {} (IDX={}), got {} (IDX={})",
                expected.1, expected.0, actual.1, actual.0,
            ),
        }
    }
}

pub(super) fn parse(s: &str) -> Result<Header, ParseError> {
    let mut builder = Header::builder();
    let mut lines = s.lines();

    let line = lines.next().ok_or(ParseError::MissingFileFormat)?;
    let file_format = parse_file_format(line)?;

    builder = builder.set_file_format(file_format);

    let mut has_header = false;

    for line in &mut lines {
        if line.starts_with("#CHROM") {
            builder = parse_header(builder, line)?;
            has_header = true;
            break;
        }

        builder = parse_record(file_format, builder, line)?;
    }

    if !has_header {
        return Err(ParseError::MissingHeader);
    }

    if lines.next().is_some() {
        return Err(ParseError::ExpectedEof);
    }

    Ok(builder.build())
}

fn parse_file_format(s: &str) -> Result<FileFormat, ParseError> {
    let record: Record = s.parse().map_err(ParseError::InvalidRecord)?;

    if record.key() != &record::key::FILE_FORMAT {
        return Err(ParseError::MissingFileFormat);
    }

    match record.value() {
        record::Value::String(value) => value.parse().map_err(ParseError::InvalidFileFormat),
        _ => Err(ParseError::InvalidRecordValue),
    }
}

fn parse_record(
    file_format: FileFormat,
    mut builder: Builder,
    line: &str,
) -> Result<Builder, ParseError> {
    use record::key;

    let record: Record = line.parse().map_err(ParseError::InvalidRecord)?;
    let (key, value) = record.into();

    builder = match key {
        key::FILE_FORMAT => {
            return Err(ParseError::UnexpectedFileFormat);
        }
        key::INFO => match value {
            record::Value::Struct(id, fields) => {
                let info = Info::try_from_fields(id, fields, file_format)
                    .map_err(ParseError::InvalidInfo)?;
                builder.add_info(info)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::FILTER => match value {
            record::Value::Struct(id, fields) => {
                let filter =
                    Filter::try_from_fields(id, fields).map_err(ParseError::InvalidFilter)?;
                builder.add_filter(filter)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::FORMAT => match value {
            record::Value::Struct(id, fields) => {
                let format = Format::try_from_fields(id, fields, file_format)
                    .map_err(ParseError::InvalidFormat)?;
                builder.add_format(format)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::ALTERNATIVE_ALLELE => match value {
            record::Value::Struct(id, fields) => {
                let alternative_allele = AlternativeAllele::try_from_fields(id, fields)
                    .map_err(ParseError::InvalidAlternativeAllele)?;
                builder.add_alternative_allele(alternative_allele)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::ASSEMBLY => match value {
            record::Value::String(value) => builder.set_assembly(value),
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::CONTIG => match value {
            record::Value::Struct(id, fields) => {
                let contig =
                    Contig::try_from_fields(id, fields).map_err(ParseError::InvalidContig)?;
                builder.add_contig(contig)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::META => match value {
            record::Value::Struct(id, fields) => {
                let meta = Meta::try_from_fields(id, fields).map_err(ParseError::InvalidMeta)?;
                builder.add_meta(meta)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::SAMPLE => match value {
            record::Value::Struct(id, fields) => {
                let sample =
                    Sample::try_from_fields(id, fields).map_err(ParseError::InvalidSample)?;
                builder.add_sample(sample)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::PEDIGREE => match value {
            record::Value::Struct(id, fields) => {
                let pedigree =
                    Pedigree::try_from_fields(id, fields).map_err(ParseError::InvalidPedigree)?;
                builder.add_pedigree(pedigree)
            }
            _ => return Err(ParseError::InvalidRecordValue),
        },
        key::PEDIGREE_DB => match value {
            record::Value::String(value) => builder.set_pedigree_db(value),
            _ => return Err(ParseError::InvalidRecordValue),
        },
        record::Key::Other(_) => {
            let record = Record::new(key, value);
            builder.insert(record)
        }
    };

    Ok(builder)
}

fn parse_header(mut builder: Builder, line: &str) -> Result<Builder, ParseError> {
    static HEADERS: &[&str] = &[
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
    ];
    static FORMAT_HEADER: &str = "FORMAT";

    let mut fields = line.split(crate::record::FIELD_DELIMITER);

    for &expected in HEADERS.iter() {
        if let Some(actual) = fields.next() {
            if actual != expected {
                return Err(ParseError::InvalidHeader(actual.into(), expected.into()));
            }
        } else {
            return Err(ParseError::InvalidHeader(String::from(""), expected.into()));
        }
    }

    if let Some(field) = fields.next() {
        if field != FORMAT_HEADER {
            return Err(ParseError::InvalidHeader(
                field.into(),
                FORMAT_HEADER.into(),
            ));
        }

        let mut sample_names = IndexSet::new();

        for sample_name in fields {
            if !sample_names.insert(sample_name.into()) {
                return Err(ParseError::DuplicateSampleName(sample_name.into()));
            }
        }

        builder = builder.set_sample_names(sample_names);
    }

    Ok(builder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_str() -> Result<(), ParseError> {
        let s = r#"##fileformat=VCFv4.3
##fileDate=20200506
##source=noodles-vcf
##assembly=file:///assemblies.fasta
##contig=<ID=sq0,length=8>
##contig=<ID=sq1,length=13>
##contig=<ID=sq2,length=21>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##FILTER=<ID=q10,Description="Quality below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=DEL,Description="Deletion">
##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
##SAMPLE=<ID=sample0,Assay=WholeGenome>
##PEDIGREE=<ID=cid,Father=fid,Mother=mid>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0
"#;

        let header = parse(s)?;

        assert_eq!(header.file_format(), FileFormat::new(4, 3));
        assert_eq!(header.infos().len(), 1);
        assert_eq!(header.filters().len(), 1);
        assert_eq!(header.formats().len(), 1);
        assert_eq!(header.alternative_alleles().len(), 1);
        assert_eq!(header.assembly(), Some("file:///assemblies.fasta"));
        assert_eq!(header.contigs().len(), 3);
        assert_eq!(header.meta().len(), 1);
        assert_eq!(header.samples().len(), 1);
        assert_eq!(header.pedigrees().len(), 1);
        assert_eq!(header.sample_names().len(), 1);

        assert_eq!(
            header.get("fileDate"),
            Some(
                &[Record::new(
                    record::Key::from("fileDate"),
                    record::Value::from("20200506"),
                )][..]
            )
        );

        assert_eq!(
            header.get("source"),
            Some(
                &[Record::new(
                    record::Key::from("source"),
                    record::Value::from("noodles-vcf"),
                )][..]
            )
        );

        Ok(())
    }

    #[test]
    fn test_from_str_without_file_format() {
        let s = r#"##ALT=<ID=DEL,Description="Deletion">
"#;

        assert_eq!(parse(s), Err(ParseError::MissingFileFormat));
    }

    #[test]
    fn test_from_str_without_assembly() -> Result<(), ParseError> {
        let s = r#"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
"#;
        let header = parse(s)?;
        assert!(header.assembly().is_none());
        Ok(())
    }

    #[test]
    fn test_from_str_with_data_after_header() {
        let s = r#"##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
##contig=<ID=sq0,length=8>
"#;

        assert_eq!(parse(s), Err(ParseError::ExpectedEof));
    }

    #[test]
    fn test_from_str_with_multiple_fileformats() {
        let s = "\
##fileformat=VCFv4.3
##fileformat=VCFv4.3
";

        assert_eq!(parse(s), Err(ParseError::UnexpectedFileFormat));
    }

    #[test]
    fn test_from_str_with_missing_headers() {
        let s = "##fileformat=VCFv4.3
";
        assert_eq!(parse(s), Err(ParseError::MissingHeader));
    }

    #[test]
    fn test_from_str_with_invalid_headers() {
        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUALITY	FILTER	INFO
";

        assert_eq!(
            parse(s),
            Err(ParseError::InvalidHeader(
                String::from("QUALITY"),
                String::from("QUAL")
            ))
        );

        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID
";

        assert_eq!(
            parse(s),
            Err(ParseError::InvalidHeader(
                String::from(""),
                String::from("REF")
            ))
        );

        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	sample0
";

        assert_eq!(
            parse(s),
            Err(ParseError::InvalidHeader(
                String::from("sample0"),
                String::from("FORMAT")
            ))
        );
    }

    #[test]
    fn test_from_str_with_duplicate_sample_names() {
        let s = "##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample0	sample0
";

        assert_eq!(
            parse(s),
            Err(ParseError::DuplicateSampleName(String::from("sample0")))
        );
    }
}
