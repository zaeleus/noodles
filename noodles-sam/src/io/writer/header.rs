mod record;

use std::io::{self, Write};

use crate::Header;
use record::{write_comment, write_program, write_read_group, write_reference_sequence};

pub(super) fn write_header<W>(writer: &mut W, header: &Header) -> io::Result<()>
where
    W: Write,
{
    if let Some(header) = header.header() {
        record::write_header(writer, header)?;
    }

    for (name, reference_sequence) in header.reference_sequences() {
        write_reference_sequence(writer, name, reference_sequence)?;
    }

    for (id, read_group) in header.read_groups() {
        write_read_group(writer, id, read_group)?;
    }

    for (id, program) in header.programs().as_ref() {
        write_program(writer, id, program)?;
    }

    for comment in header.comments() {
        write_comment(writer, comment)?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_header() -> io::Result<()> {
        use std::num::NonZeroUsize;

        const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
            Some(length) => length,
            None => unreachable!(),
        };

        const SQ1_LN: NonZeroUsize = match NonZeroUsize::new(13) {
            Some(length) => length,
            None => unreachable!(),
        };

        use crate::header::record::value::{
            Map,
            map::{self, Program, ReadGroup, ReferenceSequence, header::Version},
        };

        let header = Header::builder()
            .set_header(Map::<map::Header>::new(Version::new(1, 6)))
            .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
            .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(SQ1_LN))
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .add_read_group("rg1", Map::<ReadGroup>::default())
            .add_program("pg0", Map::<Program>::default())
            .add_program("pg1", Map::<Program>::default())
            .add_comment("noodles")
            .add_comment("sam")
            .build();

        let mut buf = Vec::new();
        write_header(&mut buf, &header)?;

        let expected = b"@HD\tVN:1.6
@SQ\tSN:sq0\tLN:8
@SQ\tSN:sq1\tLN:13
@RG\tID:rg0
@RG\tID:rg1
@PG\tID:pg0
@PG\tID:pg1
@CO\tnoodles
@CO\tsam
";

        assert_eq!(buf, expected);

        Ok(())
    }
}
