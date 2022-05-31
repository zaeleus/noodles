use std::{
    io::{self, BufRead},
    str,
};

use noodles_core::Position;

use super::read_line;
use crate::{
    alignment::Record,
    record::{Cigar, Data, Flags, MappingQuality, QualityScores, ReadName, Sequence},
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
    *record.reference_sequence_id_mut() = parse_reference_sequence_id(header, field)?;

    let field = next_field(&mut src);
    *record.alignment_start_mut() = parse_alignment_start(field)?;

    let field = next_field(&mut src);
    *record.mapping_quality_mut() = parse_mapping_quality(field)?;

    let field = next_field(&mut src);
    *record.cigar_mut() = parse_cigar(field)?;

    let field = next_field(&mut src);
    *record.mate_reference_sequence_id_mut() = parse_reference_sequence_id(header, field)?;

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

fn parse_read_name(src: &[u8]) -> io::Result<Option<ReadName>> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(None),
        _ => ReadName::try_new(src)
            .map(Some)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
    }
}

fn parse_flags(src: &[u8]) -> io::Result<Flags> {
    str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse::<u16>()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
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

fn parse_alignment_start(src: &[u8]) -> io::Result<Option<Position>> {
    str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .map(Position::new)
}

fn parse_mapping_quality(src: &[u8]) -> io::Result<Option<MappingQuality>> {
    str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
        .map(MappingQuality::new)
}

fn parse_cigar(src: &[u8]) -> io::Result<Cigar> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(Cigar::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}

fn parse_template_length(src: &[u8]) -> io::Result<i32> {
    str::from_utf8(src)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        .and_then(|s| {
            s.parse()
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
}

fn parse_sequence(src: &[u8]) -> io::Result<Sequence> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(Sequence::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}

fn parse_quality_scores(src: &[u8]) -> io::Result<QualityScores> {
    const MISSING: &[u8] = b"*";

    match src {
        MISSING => Ok(QualityScores::default()),
        _ => str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }),
    }
}

fn parse_data(src: &[u8]) -> io::Result<Data> {
    if src.is_empty() {
        Ok(Data::default())
    } else {
        str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|s| {
                s.parse()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })
    }
}
