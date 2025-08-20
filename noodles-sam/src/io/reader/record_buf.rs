pub(crate) mod cigar;
pub(crate) mod data;
mod flags;
mod mapping_quality;
mod name;
mod position;
mod quality_scores;
mod reference_sequence_id;
mod sequence;
mod template_length;

use std::{
    error, fmt,
    io::{self, BufRead},
};

use self::{
    cigar::parse_cigar, data::parse_data, mapping_quality::parse_mapping_quality, name::parse_name,
    position::parse_alignment_start, quality_scores::parse_quality_scores,
    reference_sequence_id::parse_reference_sequence_id, sequence::parse_sequence,
};
pub(crate) use self::{flags::parse_flags, template_length::parse_template_length};
use super::read_line;
use crate::{Header, alignment::RecordBuf};

pub fn read_record_buf<R>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    header: &Header,
    record: &mut RecordBuf,
) -> io::Result<usize>
where
    R: BufRead,
{
    buf.clear();

    match read_line(reader, buf)? {
        0 => Ok(0),
        n => {
            parse_record_buf(buf, header, record)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            Ok(n)
        }
    }
}

/// An error when a raw SAM record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The name is invalid.
    InvalidName(name::ParseError),
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
            Self::InvalidName(e) => Some(e),
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
            Self::InvalidName(_) => write!(f, "invalid name"),
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

pub(crate) fn parse_record_buf(
    mut src: &[u8],
    header: &Header,
    record: &mut RecordBuf,
) -> Result<(), ParseError> {
    const MISSING: &[u8] = b"*";

    match next_field(&mut src) {
        MISSING => *record.name_mut() = None,
        field => {
            parse_name(field, record.name_mut())
                .map(Some)
                .map_err(ParseError::InvalidName)?;
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

    record.cigar_mut().as_mut().clear();
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

    record.sequence_mut().as_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_sequence(field, record.sequence_mut()).map_err(ParseError::InvalidSequence)?;
    }

    record.quality_scores_mut().as_mut().clear();
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
    use std::num::NonZero;

    use super::*;

    #[test]
    fn test_parse_with_data() -> Result<(), ParseError> {
        use crate::alignment::{record::data::field::Tag, record_buf::data::field::Value};

        let header = Header::default();
        let s = b"*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\tNH:i:1\tCO:Z:ndls";
        let mut record = RecordBuf::default();
        parse_record_buf(s, &header, &mut record)?;

        let expected = RecordBuf::builder()
            .set_data(
                [
                    (Tag::ALIGNMENT_HIT_COUNT, Value::from(1)),
                    (Tag::COMMENT, Value::from("ndls")),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        assert_eq!(record, expected);

        Ok(())
    }

    #[test]
    fn test_parse_mate_reference_sequence_id() {
        use crate::header::record::value::{Map, map::ReferenceSequence};

        let header = Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::new(const { NonZero::new(13).unwrap() }),
            )
            .build();

        let reference_sequence_id = Some(0);

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"="),
            Ok(reference_sequence_id)
        );

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq0"),
            Ok(Some(0))
        );

        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq1"),
            Ok(Some(1))
        );

        assert!(matches!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"*"),
            Err(ParseError::InvalidMateReferenceSequenceId(_))
        ));

        assert!(matches!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq2"),
            Err(ParseError::InvalidMateReferenceSequenceId(_))
        ));
    }
}
