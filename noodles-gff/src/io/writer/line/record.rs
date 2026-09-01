mod attributes;
mod phase;
mod position;
mod reference_sequence_name;
mod score;
mod source;
mod strand;
mod ty;

use std::{
    borrow::Cow,
    io::{self, Write},
};

use bstr::BStr;
use percent_encoding::{AsciiSet, CONTROLS};

use self::{
    attributes::write_attributes, phase::write_phase, position::write_position,
    reference_sequence_name::write_reference_sequence_name, score::write_score,
    source::write_source, strand::write_strand, ty::write_type,
};
use crate::feature::Record;

pub(crate) fn write_record<W, R>(writer: &mut W, record: &R) -> io::Result<()>
where
    W: Write,
    R: Record + ?Sized,
{
    write_reference_sequence_name(writer, record.reference_sequence_name().as_ref())?;

    write_separator(writer)?;
    write_source(writer, record.source().as_ref())?;

    write_separator(writer)?;
    let ty = record.ty();
    write_type(writer, ty.as_ref())?;

    write_separator(writer)?;
    write_position(writer, record.feature_start()?)?;

    write_separator(writer)?;
    write_position(writer, record.feature_end()?)?;

    write_separator(writer)?;
    write_score(writer, record.score().transpose()?)?;

    write_separator(writer)?;
    write_strand(writer, record.strand()?)?;

    write_separator(writer)?;
    write_phase(writer, ty.as_ref(), record.phase().transpose()?)?;

    write_separator(writer)?;
    write_attributes(writer, record.attributes().as_ref())?;

    Ok(())
}

fn write_missing<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const MISSING: u8 = b'.';
    writer.write_all(&[MISSING])
}

fn write_separator<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    const SEPARATOR: u8 = b'\t';
    writer.write_all(&[SEPARATOR])
}

// § "Description of the Format" (2020-08-18): "Literal use of tab, newline, carriage return, the
// percent (%) sign, and control characters must be encoded using RFC 3986 Percent-Encoding; no
// other characters may be encoded."
fn percent_encode(s: &BStr) -> Cow<'_, str> {
    const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS.add(b'%');
    percent_encoding::percent_encode(s, PERCENT_ENCODE_SET).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::feature::RecordBuf;

    #[test]
    fn test_write_record() -> io::Result<()> {
        let mut buf = Vec::new();
        let record = RecordBuf::default();
        write_record(&mut buf, &record)?;
        assert_eq!(buf, b".\t.\t.\t1\t1\t.\t.\t.\t.");
        Ok(())
    }

    #[test]
    fn test_percent_encode() {
        assert_eq!(percent_encode(BStr::new("")), "");
        assert_eq!(percent_encode(BStr::new("noodles")), "noodles");
        assert_eq!(percent_encode(BStr::new("\t\n\r%\0")), "%09%0A%0D%25%00");
    }
}
