mod cigar;
pub(crate) mod data;
mod quality_scores;
mod reference_sequence_id;
mod sequence;

pub(crate) use self::{
    cigar::parse_cigar, data::parse_data, quality_scores::parse_quality_scores,
    sequence::parse_sequence,
};

use std::{
    error, fmt,
    io::{self, BufRead},
};

use noodles_core::Position;

use self::reference_sequence_id::parse_reference_sequence_id;
use super::read_line;
use crate::{
    alignment::Record,
    record::{Flags, MappingQuality, ReadName},
    Header,
};

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
            parse_record(buf, header, record)?;
            Ok(n)
        }
    }
}

/// An error when a raw SAM record fails to parse.
#[allow(clippy::enum_variant_names)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseError {
    /// The reference sequence ID is invalid.
    InvalidReferenceSequenceId(reference_sequence_id::ParseError),
    /// The CIGAR is invalid.
    InvalidCigar(cigar::ParseError),
    /// The mate reference sequence ID is invalid.
    InvalidMateReferenceSequenceId(reference_sequence_id::ParseError),
    /// The sequence is invalid.
    InvalidSequence(sequence::ParseError),
    /// The data is invalid.
    InvalidData(data::ParseError),
}

impl error::Error for ParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidReferenceSequenceId(e) => Some(e),
            Self::InvalidCigar(e) => Some(e),
            Self::InvalidMateReferenceSequenceId(e) => Some(e),
            Self::InvalidSequence(e) => Some(e),
            Self::InvalidData(e) => Some(e),
        }
    }
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidReferenceSequenceId(_) => write!(f, "invalid reference sequence ID"),
            Self::InvalidCigar(_) => write!(f, "invalid CIGAR"),
            Self::InvalidMateReferenceSequenceId(_) => {
                write!(f, "invalid mate reference sequence ID")
            }
            Self::InvalidSequence(_) => write!(f, "invalid sequence"),
            Self::InvalidData(_) => write!(f, "invalid data"),
        }
    }
}

pub(crate) fn parse_record(mut src: &[u8], header: &Header, record: &mut Record) -> io::Result<()> {
    const MISSING: &[u8] = b"*";

    *record.read_name_mut() = match next_field(&mut src) {
        MISSING => None,
        field => parse_read_name(field).map(Some)?,
    };

    let field = next_field(&mut src);
    *record.flags_mut() = parse_flags(field)?;

    let reference_sequence_id = match next_field(&mut src) {
        MISSING => None,
        field => parse_reference_sequence_id(header, field)
            .map(Some)
            .map_err(ParseError::InvalidReferenceSequenceId)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    };

    *record.reference_sequence_id_mut() = reference_sequence_id;

    let field = next_field(&mut src);
    *record.alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.mapping_quality_mut() = parse_mapping_quality(field)?;

    record.cigar_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_cigar(field, record.cigar_mut())
            .map_err(ParseError::InvalidCigar)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    *record.mate_reference_sequence_id_mut() = match next_field(&mut src) {
        MISSING => None,
        field => parse_mate_reference_sequence_id(header, reference_sequence_id, field)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    };

    let field = next_field(&mut src);
    *record.mate_alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.template_length_mut() = parse_template_length(field)?;

    record.sequence_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_sequence(field, record.sequence_mut())
            .map_err(ParseError::InvalidSequence)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    }

    record.quality_scores_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_quality_scores(field, record.quality_scores_mut())?;
    }

    record.data_mut().clear();
    let field = next_field(&mut src);
    parse_data(field, record.data_mut())
        .map_err(ParseError::InvalidData)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

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

pub(crate) fn parse_read_name(src: &[u8]) -> io::Result<ReadName> {
    ReadName::try_new(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

pub(crate) fn parse_flags(src: &[u8]) -> io::Result<Flags> {
    lexical_core::parse::<u16>(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Flags::from)
}

pub(crate) fn parse_alignment_start(src: &[u8]) -> io::Result<Option<Position>> {
    lexical_core::parse(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Position::new)
}

pub(crate) fn parse_mapping_quality(src: &[u8]) -> io::Result<Option<MappingQuality>> {
    lexical_core::parse(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(MappingQuality::new)
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

pub(crate) fn parse_template_length(src: &[u8]) -> io::Result<i32> {
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use super::*;

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
