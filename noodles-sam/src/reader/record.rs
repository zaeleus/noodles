pub(crate) mod cigar;
pub(crate) mod data;
mod flags;
mod mapping_quality;
mod position;
mod quality_scores;
mod read_name;
mod reference_sequence_id;
mod sequence;
mod template_length;

pub(crate) use self::{
    cigar::parse_cigar, flags::parse_flags, mapping_quality::parse_mapping_quality,
    position::parse_alignment_start, read_name::parse_read_name, sequence::parse_sequence,
    template_length::parse_template_length,
};

use std::{
    error, fmt,
    io::{self, BufRead},
};

use self::{
    data::parse_data, quality_scores::parse_quality_scores,
    reference_sequence_id::parse_reference_sequence_id,
};
use super::read_line;
use crate::{alignment::Record, Header};

pub fn read_record<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    header: &Header,
    record: &mut Record,
) -> io::Result<usize>
where
    R: BufRead,
{
    buf.clear();

    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            parse_record(buf, header, record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}

/// An error when a raw SAM record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The read name is invalid.
    InvalidReadName(read_name::ParseError),
    /// The flags are invalid.
    InvalidFlags(flags::ParseError),
    /// The reference sequence ID is invalid.
    InvalidReferenceSequenceId(reference_sequence_id::ParseError),
    /// The position is invalid.
    InvalidPosition(position::ParseError),
    /// The mapping quality is invalid.
    InvalidMappingQuality(mapping_quality::ParseError),
    /// The CIGAR is invalid.
    InvalidCigar(cigar::ParseError),
    /// The mate reference sequence ID is invalid.
    InvalidMateReferenceSequenceId(reference_sequence_id::ParseError),
    /// The mate position is invalid.
    InvalidMatePosition(position::ParseError),
    /// The template length is invalid.
    InvalidTemplateLength(template_length::ParseError),
    /// The sequence is invalid.
    InvalidSequence(sequence::ParseError),
    /// The quality scores are invalid.
    InvalidQualityScores(quality_scores::ParseError),
    /// The data is invalid.
    InvalidData(data::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidReadName(e) => Some(e),
            Self::InvalidFlags(e) => Some(e),
            Self::InvalidReferenceSequenceId(e) => Some(e),
            Self::InvalidPosition(e) => Some(e),
            Self::InvalidMappingQuality(e) => Some(e),
            Self::InvalidCigar(e) => Some(e),
            Self::InvalidMateReferenceSequenceId(e) => Some(e),
            Self::InvalidMatePosition(e) => Some(e),
            Self::InvalidTemplateLength(e) => Some(e),
            Self::InvalidSequence(e) => Some(e),
            Self::InvalidQualityScores(e) => Some(e),
            Self::InvalidData(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidReadName(_) => write!(f, "invalid read name"),
            Self::InvalidFlags(_) => write!(f, "invalid flags"),
            Self::InvalidReferenceSequenceId(_) => write!(f, "invalid reference sequence ID"),
            Self::InvalidPosition(_) => write!(f, "invalid position"),
            Self::InvalidMappingQuality(_) => write!(f, "invalid mapping quality"),
            Self::InvalidCigar(_) => write!(f, "invalid CIGAR"),
            Self::InvalidMateReferenceSequenceId(_) => {
                write!(f, "invalid mate reference sequence ID")
            }
            Self::InvalidMatePosition(_) => write!(f, "invalid mate position"),
            Self::InvalidTemplateLength(_) => write!(f, "invalid template length"),
            Self::InvalidSequence(_) => write!(f, "invalid sequence"),
            Self::InvalidQualityScores(_) => write!(f, "invalid quality scores"),
            Self::InvalidData(_) => write!(f, "invalid data"),
        }
    }
}

pub(crate) fn parse_record(
    mut src: &[u8],
    header: &Header,
    record: &mut Record,
) -> Result<(), ParseError> {
    const MISSING: &[u8] = b"*";

    match next_field(&mut src) {
        MISSING => *record.read_name_mut() = None,
        field => {
            parse_read_name(field, record.read_name_mut())
                .map(Some)
                .map_err(ParseError::InvalidReadName)?;
        }
    };

    let field = next_field(&mut src);
    *record.flags_mut() = parse_flags(field).map_err(ParseError::InvalidFlags)?;

    let reference_sequence_id = match next_field(&mut src) {
        MISSING => None,
        field => parse_reference_sequence_id(header, field)
            .map(Some)
            .map_err(ParseError::InvalidReferenceSequenceId)?,
    };

    *record.reference_sequence_id_mut() = reference_sequence_id;

    let field = next_field(&mut src);
    *record.alignment_start_mut() =
        parse_alignment_start(field).map_err(ParseError::InvalidPosition)?;

    let field = next_field(&mut src);
    *record.mapping_quality_mut() =
        parse_mapping_quality(field).map_err(ParseError::InvalidMappingQuality)?;

    record.cigar_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_cigar(field, record.cigar_mut()).map_err(ParseError::InvalidCigar)?;
    }

    *record.mate_reference_sequence_id_mut() = match next_field(&mut src) {
        MISSING => None,
        field => parse_mate_reference_sequence_id(header, reference_sequence_id, field)?,
    };

    let field = next_field(&mut src);
    *record.mate_alignment_start_mut() =
        parse_alignment_start(field).map_err(ParseError::InvalidMatePosition)?;

    let field = next_field(&mut src);
    *record.template_length_mut() =
        parse_template_length(field).map_err(ParseError::InvalidTemplateLength)?;

    record.sequence_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_sequence(field, record.sequence_mut()).map_err(ParseError::InvalidSequence)?;
    }

    record.quality_scores_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_quality_scores(field, record.sequence().len(), record.quality_scores_mut())
            .map_err(ParseError::InvalidQualityScores)?;
    }

    record.data_mut().clear();
    parse_data(src, record.data_mut()).map_err(ParseError::InvalidData)?;

    Ok(())
}

fn next_field<'a>(src: &mut &'a [u8]) -> &'a [u8] {
    use memchr::memchr;

    const DELIMITER: u8 = b'\t';

    let (field, rest) = if let Some(i) = memchr(DELIMITER, src) {
        let (field, rest) = src.split_at(i);
        (field, &rest[1..])
    } else {
        src.split_at(src.len())
    };

    *src = rest;

    field
}

fn parse_mate_reference_sequence_id(
    header: &Header,
    reference_sequence_id: Option<usize>,
    src: &[u8],
) -> Result<Option<usize>, ParseError> {
    const EQ: &[u8] = b"=";

    match src {
        EQ => Ok(reference_sequence_id),
        _ => parse_reference_sequence_id(header, src)
            .map(Some)
            .map_err(ParseError::InvalidMateReferenceSequenceId),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_with_data() -> Result<(), ParseError> {
        use crate::record::data::field::{tag, Value};

        let header = Header::default();
        let s = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tNH:i:1\tCO:Z:ndls";
        let mut record = Record::default();
        parse_record(s, &header, &mut record)?;

        let expected = Record::builder()
            .set_data(
                [
                    (tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
                    (tag::COMMENT, Value::String(String::from("ndls"))),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        assert_eq!(record, expected);

        Ok(())
    }

    #[test]
    fn test_parse_mate_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZeroUsize;

        use crate::header::record::value::{map::ReferenceSequence, Map};

        let header = Header::builder()
            .add_reference_sequence(
                "sq0".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
            )
            .add_reference_sequence(
                "sq1".parse()?,
                Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
            )
            .build();

        let reference_sequence_id = Some(0);

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"=")?,
            reference_sequence_id
        );

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq0")?,
            Some(0)
        );

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq1")?,
            Some(1)
        );

        assert!(matches!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"*"),
            Err(ParseError::InvalidMateReferenceSequenceId(_))
        ));

        assert!(matches!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq2"),
            Err(ParseError::InvalidMateReferenceSequenceId(_))
        ));

        Ok(())
    }
}
