use std::{
    error, fmt,
    io::{self, BufRead, Lines},
};

use noodles_core::{position, region::Interval, Position, Region};

use crate::index::{header::format::CoordinateSystem, Header};

/// An iterator that filters raw records that intersect the given region.
pub struct FilteredLines<'h, R> {
    lines: Lines<R>,
    header: &'h Header,
    region: Region,
}

impl<'h, R> FilteredLines<'h, R>
where
    R: BufRead,
{
    /// Creates a filtered lines iterator.
    pub fn new(reader: R, header: &'h Header, region: Region) -> Self {
        Self {
            lines: reader.lines(),
            header,
            region,
        }
    }
}

impl<'h, R> Iterator for FilteredLines<'h, R>
where
    R: BufRead,
{
    type Item = io::Result<String>;

    fn next(&mut self) -> Option<Self::Item> {
        let line_comment_prefix = char::from(self.header.line_comment_prefix());

        loop {
            let line = match self.lines.next()? {
                Ok(s) => s,
                Err(e) => return Some(Err(e)),
            };

            if line.starts_with(line_comment_prefix) {
                continue;
            }

            let (reference_sequence_name, interval) = match parse_record(&line, self.header) {
                Ok(record) => record,
                Err(e) => return Some(Err(io::Error::new(io::ErrorKind::InvalidData, e))),
            };

            if intersects(reference_sequence_name, interval, &self.region) {
                return Some(Ok(line));
            }
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParseRecordError {
    /// The reference sequence name is missing.
    MissingReferenceSequenceName,
    /// The start position is missing.
    MissingStartPosition,
    /// The start position is invalid.
    InvalidStartPosition(ParsePositionError),
    /// The end position is missing.
    MissingEndPosition,
    /// The end position is invalid.
    InvalidEndPosition(ParsePositionError),
}

impl error::Error for ParseRecordError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::InvalidStartPosition(e) => Some(e),
            Self::InvalidEndPosition(e) => Some(e),
            _ => None,
        }
    }
}

impl fmt::Display for ParseRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingReferenceSequenceName => write!(f, "missing reference sequence name"),
            Self::MissingStartPosition => write!(f, "missing start position"),
            Self::InvalidStartPosition(_) => write!(f, "invalid start position"),
            Self::MissingEndPosition => write!(f, "missing end position"),
            Self::InvalidEndPosition(_) => write!(f, "invalid end position"),
        }
    }
}

fn parse_record<'a>(s: &'a str, header: &Header) -> Result<(&'a str, Interval), ParseRecordError> {
    const DELIMITER: char = '\t';

    let fields: Vec<_> = s.split(DELIMITER).collect();

    let reference_sequence_name = fields
        .get(header.reference_sequence_name_index())
        .ok_or(ParseRecordError::MissingReferenceSequenceName)?;

    let raw_start = fields
        .get(header.start_position_index())
        .ok_or(ParseRecordError::MissingStartPosition)?;

    let coordinate_system = header.format().coordinate_system();
    let start = parse_start_position(raw_start, coordinate_system)
        .map_err(ParseRecordError::InvalidStartPosition)?;

    let end = if let Some(i) = header.end_position_index() {
        fields
            .get(i)
            .ok_or(ParseRecordError::MissingEndPosition)
            .and_then(|s| {
                s.parse()
                    .map_err(ParsePositionError::Parse)
                    .map_err(ParseRecordError::InvalidEndPosition)
            })?
    } else {
        // _The Tabix index file format_: "Field `col_beg` may equal `col_end`, and in this case,
        // the end of a region is `end=beg+1`."
        start.checked_add(1).expect("attempt to add with overflow")
    };

    Ok((reference_sequence_name, Interval::from(start..=end)))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ParsePositionError {
    /// The input is invalid.
    Parse(position::ParseError),
    /// The position is invalid.
    Invalid(position::TryFromIntError),
}

impl error::Error for ParsePositionError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Parse(e) => Some(e),
            Self::Invalid(e) => Some(e),
        }
    }
}

impl fmt::Display for ParsePositionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(_) => write!(f, "invalid input"),
            Self::Invalid(_) => write!(f, "invalid input"),
        }
    }
}

fn parse_start_position(
    s: &str,
    coordinate_system: CoordinateSystem,
) -> Result<Position, ParsePositionError> {
    match coordinate_system {
        CoordinateSystem::Gff => s.parse().map_err(ParsePositionError::Parse),
        CoordinateSystem::Bed => s
            .parse::<usize>()
            .map_err(ParsePositionError::Parse)
            .and_then(|n| Position::try_from(n + 1).map_err(ParsePositionError::Invalid)),
    }
}

fn intersects(reference_sequence_name: &str, interval: Interval, region: &Region) -> bool {
    reference_sequence_name == region.name() && interval.intersects(region.interval())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_next() -> Result<(), Box<dyn std::error::Error>> {
        let data = b"sq0\t8\t13
sq1\t8\t13
sq1\t21\t34
# noodles
sq1\t34\t55
sq1\t89\t144
sq2\t8\t13
";

        let reader = &data[..];

        let header = crate::index::Header::builder()
            .set_start_position_index(1)
            .set_end_position_index(Some(2))
            .set_reference_sequence_names(
                [
                    String::from("sq0"),
                    String::from("sq1"),
                    String::from("sq2"),
                ]
                .into_iter()
                .collect(),
            )
            .build();

        let region = "sq1:21-55".parse()?;

        let actual: Vec<_> =
            FilteredLines::new(reader, &header, region).collect::<Result<_, _>>()?;
        let expected = vec!["sq1\t21\t34", "sq1\t34\t55"];
        assert_eq!(actual, expected);

        Ok(())
    }
}
