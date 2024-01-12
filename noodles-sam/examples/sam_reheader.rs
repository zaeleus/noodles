//! Replaces the header of a SAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{env, io};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = sam::io::reader::Builder::default().build_from_path(src)?;
    let mut header = reader.read_header()?;

    header.add_comment("a comment added by noodles-sam");

    let stdout = io::stdout().lock();
    let mut writer = sam::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
