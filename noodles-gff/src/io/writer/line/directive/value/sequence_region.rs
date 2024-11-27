use std::io::{self, Write};

use crate::directive_buf::value::SequenceRegion;

pub(crate) fn write_sequence_region<W>(
    writer: &mut W,
    sequence_region: &SequenceRegion,
) -> io::Result<()>
where
    W: Write,
{
    write!(
        writer,
        "{} {} {}",
        sequence_region.reference_sequence_name(),
        sequence_region.start(),
        sequence_region.end()
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_sequence_region() -> io::Result<()> {
        let mut buf = Vec::new();
        let sequence_region = SequenceRegion::new(String::from("sq0"), 8, 13);
        write_sequence_region(&mut buf, &sequence_region)?;
        assert_eq!(buf, b"sq0 8 13");
        Ok(())
    }
}
