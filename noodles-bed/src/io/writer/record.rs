mod feature_end;
mod feature_start;
mod name;
mod other_fields;
mod reference_sequence_name;
mod score;
mod strand;

use std::io::{self, Write};

use self::{
    feature_end::write_feature_end, feature_start::write_feature_start, name::write_name,
    other_fields::write_other_fields, reference_sequence_name::write_reference_sequence_name,
    score::write_score, strand::write_strand,
};
use crate::feature::Record;

pub(super) fn write_record_3<W, R, const N: usize>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record<N>,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    let feature_start = record.feature_start()?;
    write_feature_start(writer, feature_start)?;

    write_separator(writer)?;
    let feature_end = record.feature_end().transpose()?;
    write_feature_end(writer, feature_end)?;

    write_other_fields(writer, record.other_fields().as_ref())?;

    write_newline(writer)?;

    Ok(())
}

pub(super) fn write_record_4<W, R, const N: usize>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record<N>,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    let feature_start = record.feature_start()?;
    write_feature_start(writer, feature_start)?;

    write_separator(writer)?;
    let feature_end = record.feature_end().transpose()?;
    write_feature_end(writer, feature_end)?;

    write_separator(writer)?;
    let name = record
        .name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing name"))?;
    write_name(writer, name)?;

    write_other_fields(writer, record.other_fields().as_ref())?;

    write_newline(writer)?;

    Ok(())
}

pub(super) fn write_record_5<W, R, const N: usize>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record<N>,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    let feature_start = record.feature_start()?;
    write_feature_start(writer, feature_start)?;

    write_separator(writer)?;
    let feature_end = record.feature_end().transpose()?;
    write_feature_end(writer, feature_end)?;

    write_separator(writer)?;
    let name = record
        .name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing name"))?;
    write_name(writer, name)?;

    write_separator(writer)?;
    let score = record
        .score()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing score"))?;
    write_score(writer, score)?;

    write_other_fields(writer, record.other_fields().as_ref())?;

    write_newline(writer)?;

    Ok(())
}

pub(super) fn write_record_6<W, R, const N: usize>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record<N>,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    let feature_start = record.feature_start()?;
    write_feature_start(writer, feature_start)?;

    write_separator(writer)?;
    let feature_end = record.feature_end().transpose()?;
    write_feature_end(writer, feature_end)?;

    write_separator(writer)?;
    let name = record
        .name()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing name"))?;
    write_name(writer, name)?;

    write_separator(writer)?;
    let score = record
        .score()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing score"))?;
    write_score(writer, score)?;

    write_separator(writer)?;
    let strand = record
        .strand()
        .transpose()?
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "missing score"))?;
    write_strand(writer, strand)?;

    write_other_fields(writer, record.other_fields().as_ref())?;

    write_newline(writer)?;

    Ok(())
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'\t';
    writer.write_all(&[SEPARATOR])
}

fn write_newline<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const LINE_FEED: u8 = b'\n';
    writer.write_all(&[LINE_FEED])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_record_3() -> io::Result<()> {
        let mut buf = Vec::new();
        let record = crate::Record::<3>::default();
        write_record_3(&mut buf, &record)?;
        assert_eq!(buf, b"sq0\t0\t1\n");
        Ok(())
    }
}
