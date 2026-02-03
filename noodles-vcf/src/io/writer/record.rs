mod alternate_bases;
mod filters;
mod ids;
mod info;
mod position;
mod quality_score;
mod reference_bases;
mod reference_sequence_name;
mod samples;

use std::{
    error, fmt,
    io::{self, Write},
};

use self::{
    alternate_bases::write_alternate_bases, filters::write_filters, ids::write_ids,
    info::write_info, position::write_position, quality_score::write_quality_score,
    reference_bases::write_reference_bases, reference_sequence_name::write_reference_sequence_name,
    samples::write_samples,
};
use crate::{Header, variant::Record};

const MISSING: &[u8] = b".";

/// An error returns when a record fails to write.
#[derive(Debug)]
pub enum WriteError {
    // I/O error.
    Io(io::Error),
    // The reference sequence name is invalid.
    InvalidReferenceSequenceName(reference_sequence_name::WriteError),
    // The IDs are invalid.
    InvalidIds(ids::WriteError),
    // The reference bases are invalid.
    InvalidReferenceBases(reference_bases::WriteError),
    // The alternate bases are invalid.
    InvalidAlternateBases(alternate_bases::WriteError),
    // The filters are invalid.
    InvalidFilters(filters::WriteError),
    // The info fields are invalid.
    InvalidInfo(info::WriteError),
    // The samples are invalid.
    InvalidSamples(samples::WriteError),
}

impl error::Error for WriteError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            Self::InvalidReferenceSequenceName(e) => Some(e),
            Self::InvalidIds(e) => Some(e),
            Self::InvalidReferenceBases(e) => Some(e),
            Self::InvalidAlternateBases(e) => Some(e),
            Self::InvalidFilters(e) => Some(e),
            Self::InvalidInfo(e) => Some(e),
            Self::InvalidSamples(e) => Some(e),
        }
    }
}

impl fmt::Display for WriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(_) => write!(f, "I/O error"),
            Self::InvalidReferenceSequenceName(_) => write!(f, "invalid reference sequence name"),
            Self::InvalidIds(_) => write!(f, "invalid IDs"),
            Self::InvalidReferenceBases(_) => write!(f, "invalid reference bases"),
            Self::InvalidAlternateBases(_) => write!(f, "invalid alternate bases"),
            Self::InvalidFilters(_) => write!(f, "invalid filters"),
            Self::InvalidInfo(_) => write!(f, "invalid info"),
            Self::InvalidSamples(_) => write!(f, "invalid samples"),
        }
    }
}

pub(super) fn write_record<W, R>(
    writer: &mut W,
    header: &Header,
    record: &R,
) -> Result<(), WriteError>
where
    W: Write,
    R: Record + ?Sized,
{
    let reference_sequence_name = record
        .reference_sequence_name(header)
        .map_err(WriteError::Io)?;

    write_reference_sequence_name(writer, reference_sequence_name)
        .map_err(WriteError::InvalidReferenceSequenceName)?;

    write_separator(writer)?;
    let position = record.variant_start().transpose().map_err(WriteError::Io)?;
    write_position(writer, position).map_err(WriteError::Io)?;

    write_separator(writer)?;
    write_ids(writer, record.ids()).map_err(WriteError::InvalidIds)?;

    write_separator(writer)?;
    write_reference_bases(writer, record.reference_bases())
        .map_err(WriteError::InvalidReferenceBases)?;

    write_separator(writer)?;
    write_alternate_bases(writer, record.alternate_bases())
        .map_err(WriteError::InvalidAlternateBases)?;

    write_separator(writer)?;
    let quality_score = record.quality_score().transpose().map_err(WriteError::Io)?;
    write_quality_score(writer, quality_score).map_err(WriteError::Io)?;

    write_separator(writer)?;
    write_filters(writer, header, record.filters()).map_err(WriteError::InvalidFilters)?;

    write_separator(writer)?;
    write_info(writer, header, record.info()).map_err(WriteError::InvalidInfo)?;

    let samples = record.samples().map_err(WriteError::Io)?;

    if !samples.is_empty() {
        write_separator(writer)?;
        write_samples(writer, header, samples).map_err(WriteError::InvalidSamples)?;
    }

    writer.write_all(b"\n").map_err(WriteError::Io)?;

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> Result<(), WriteError>
where
    W: Write,
{
    const SEPARATOR: &[u8] = b"\t";
    writer.write_all(SEPARATOR).map_err(WriteError::Io)
}

#[cfg(test)]
mod tests {
    use noodles_core::Position;

    use super::*;
    use crate::variant::RecordBuf;

    #[test]
    fn test_write_record() -> Result<(), WriteError> {
        let record = RecordBuf::builder()
            .set_reference_sequence_name("sq0")
            .set_variant_start(Position::MIN)
            .set_reference_bases("A")
            .build();

        let header = Header::default();
        let mut buf = Vec::new();
        write_record(&mut buf, &header, &record)?;
        assert_eq!(buf, b"sq0\t1\t.\tA\t.\t.\t.\t.\n");

        Ok(())
    }
}
