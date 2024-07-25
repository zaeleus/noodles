#![allow(dead_code)]

mod feature_end;
mod feature_start;
mod other_fields;
mod reference_sequence_name;
mod score;
mod strand;

use std::io::{self, Write};

use self::{
    feature_end::write_feature_end, feature_start::write_feature_start,
    other_fields::write_other_fields, reference_sequence_name::write_reference_sequence_name,
};
use crate::Record;

pub(super) fn write_record<W, const N: usize>(writer: &mut W, record: &Record<N>) -> io::Result<()>
where
    W: Write,
{
    write_reference_sequence_name(writer, record.reference_sequence_name())?;

    write_separator(writer)?;
    let feature_start = record.feature_start()?;
    write_feature_start(writer, feature_start)?;

    write_separator(writer)?;
    let feature_end = record.feature_end().transpose()?;
    write_feature_end(writer, feature_end)?;

    write_other_fields(writer, &record.other_fields())?;

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
    fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();
        let record = Record::<3>::default();
        write_record(&mut buf, &record)?;
        assert_eq!(buf, b"sq0\t0\t1\n");
        Ok(())
    }
}
