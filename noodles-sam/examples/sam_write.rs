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

const SQ0_LN: NonZeroUsize = match NonZeroUsize::new(8) {
    Some(length) => length,
    None => unreachable!(),
};

const SQ1_LN: NonZeroUsize = match NonZeroUsize::new(13) {
    Some(length) => length,
    None => unreachable!(),
};

const SQ2_LN: NonZeroUsize = match NonZeroUsize::new(21) {
    Some(length) => length,
    None => unreachable!(),
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    let header = sam::Header::builder()
        .set_header(Default::default())
        .add_reference_sequence("sq0", Map::<ReferenceSequence>::new(SQ0_LN))
        .add_reference_sequence("sq1", Map::<ReferenceSequence>::new(SQ1_LN))
        .add_reference_sequence("sq2", Map::<ReferenceSequence>::new(SQ2_LN))
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    writer.write_header(&header)?;

    let record = RecordBuf::default();
    writer.write_alignment_record(&header, &record)?;

    Ok(())
}
