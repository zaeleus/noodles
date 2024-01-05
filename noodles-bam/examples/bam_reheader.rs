//! Replaces the SAM header of a BAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{env, io};

use noodles_bam as bam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::reader::Builder.build_from_path(src)?;

    let mut header = reader.read_header()?;
    header.add_comment("a comment added by noodles-bam");

    let stdout = io::stdout().lock();
    let mut writer = bam::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.record_bufs(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
