//! Creates a new SAM file.
//!
//! This writes a SAM header and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{io, num::NonZeroUsize};

use noodles_sam::{
    self as sam,
    alignment::{io::Write, RecordBuf},
    header::record::value::{
        map::{Program, ReferenceSequence},
        Map,
    },
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    let header = sam::Header::builder()
        .set_header(Default::default())
        .add_reference_sequence(
            "sq0".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
        )
        .add_reference_sequence(
            "sq1".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
        )
        .add_reference_sequence(
            "sq2".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(21)?),
        )
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    writer.write_header(&header)?;

    let record = RecordBuf::default();
    writer.write_alignment_record(&header, &record)?;

    Ok(())
}
