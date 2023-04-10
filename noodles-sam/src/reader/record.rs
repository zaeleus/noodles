mod cigar;
pub(crate) mod data;
mod quality_scores;
mod reference_sequence_id;
mod sequence;

pub(crate) use self::{
    cigar::parse_cigar, data::parse_data, quality_scores::parse_quality_scores,
    sequence::parse_sequence,
};

use std::io::{self, BufRead};

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
        field => parse_reference_sequence_id(header, field).map(Some)?,
    };

    *record.reference_sequence_id_mut() = reference_sequence_id;

    let field = next_field(&mut src);
    *record.alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.mapping_quality_mut() = parse_mapping_quality(field)?;

    record.cigar_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_cigar(field, record.cigar_mut())?;
    }

    *record.mate_reference_sequence_id_mut() = match next_field(&mut src) {
        MISSING => None,
        field => parse_mate_reference_sequence_id(header, reference_sequence_id, field)?,
    };

    let field = next_field(&mut src);
    *record.mate_alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.template_length_mut() = parse_template_length(field)?;

    record.sequence_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_sequence(field, record.sequence_mut())?;
    }

    record.quality_scores_mut().clear();
    let field = next_field(&mut src);
    if field != MISSING {
        parse_quality_scores(field, record.quality_scores_mut())?;
    }

    record.data_mut().clear();
    let field = next_field(&mut src);
    parse_data(field, record.data_mut())?;

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
) -> io::Result<Option<usize>> {
    const EQ: &[u8] = b"=";

    match src {
        EQ => Ok(reference_sequence_id),
        _ => parse_reference_sequence_id(header, src).map(Some),
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
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        assert!(matches!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, b"sq2"),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
