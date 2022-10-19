//! Replaces the header of a SAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    let mut header: sam::Header = reader.read_header()?.parse()?;

    header.add_comment("a comment added by noodles-sam");

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
