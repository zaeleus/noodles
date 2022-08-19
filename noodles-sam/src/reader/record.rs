mod cigar;
pub(crate) mod data;
mod quality_scores;
mod sequence;

pub(crate) use self::{
    cigar::parse_cigar, data::parse_data, quality_scores::parse_quality_scores,
    sequence::parse_sequence,
};

use std::{
    io::{self, BufRead},
    str,
};

use noodles_core::Position;

use super::read_line;
use crate::{
    alignment::Record,
    record::{Flags, MappingQuality, ReadName},
    Header,
};

pub fn read_record<R>(reader: &mut R, header: &Header, record: &mut Record) -> io::Result<usize>
where
    R: BufRead,
{
    let mut buf = Vec::new();

    match read_line(reader, &mut buf)? {
        0 => Ok(0),
        n => {
            parse_record(&buf, header, record)?;
            Ok(n)
        }
    }
}

pub(crate) fn parse_record(mut src: &[u8], header: &Header, record: &mut Record) -> io::Result<()> {
    let field = next_field(&mut src);
    *record.read_name_mut() = parse_read_name(field)?;

    let field = next_field(&mut src);
    *record.flags_mut() = parse_flags(field)?;

    let field = next_field(&mut src);
    let reference_sequence_id = parse_reference_sequence_id(header, field)?;
    *record.reference_sequence_id_mut() = reference_sequence_id;

    let field = next_field(&mut src);
    *record.alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.mapping_quality_mut() = parse_mapping_quality(field)?;

    let field = next_field(&mut src);
    *record.cigar_mut() = parse_cigar(field)?;

    let field = next_field(&mut src);
    *record.mate_reference_sequence_id_mut() =
        parse_mate_reference_sequence_id(header, reference_sequence_id, field)?;

    let field = next_field(&mut src);
    *record.mate_alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.template_length_mut() = parse_template_length(field)?;

    let field = next_field(&mut src);
    *record.sequence_mut() = parse_sequence(field)?;

    let field = next_field(&mut src);
    *record.quality_scores_mut() = parse_quality_scores(field)?;

    let field = next_field(&mut src);
    *record.data_mut() = parse_data(field)?;

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

pub(crate) fn parse_read_name(src: &[u8]) -> io::Result<Option<ReadName>> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(None),
        _ => ReadName::try_new(src)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

pub(crate) fn parse_flags(src: &[u8]) -> io::Result<Flags> {
    lexical_core::parse::<u16>(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .map(Flags::from)
}

fn parse_reference_sequence_id(header: &Header, src: &[u8]) -> io::Result<Option<usize>> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(None),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                header
                    .reference_sequences()
                    .get_index_of(s)
                    .map(Some)
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "invalid reference sequence name",
                        )
                    })
            }),
    }
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
        _ => parse_reference_sequence_id(header, src),
    }
}

pub(crate) fn parse_template_length(src: &[u8]) -> io::Result<i32> {
    lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

#[cfg(test)]
mod tests {
    use crate::header::record::value::{map::ReferenceSequence, Map};

    use super::*;

    #[test]
    fn test_parse_mate_reference_sequence_id() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .add_reference_sequence(Map::<ReferenceSequence>::new("sq0".parse()?, 8)?)
            .add_reference_sequence(Map::<ReferenceSequence>::new("sq1".parse()?, 13)?)
            .build();

        let reference_sequence_id = Some(0);

        let src = b"*";
        assert!(parse_mate_reference_sequence_id(&header, reference_sequence_id, src)?.is_none());

        let src = b"=";
        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, src)?,
            reference_sequence_id
        );

        let src = b"sq1";
        assert_eq!(
            parse_mate_reference_sequence_id(&header, reference_sequence_id, src)?,
            Some(1)
        );

        Ok(())
    }
}
