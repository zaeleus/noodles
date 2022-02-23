//! Prints the header of a BAM file.
//!
//! A BAM file header is a SAM header.
//!
//! The result matches the output of `samtools head <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let raw_header = reader.read_header()?;

    if raw_header.is_empty() {
        let reference_sequences = reader.read_reference_sequences()?;

        let header = sam::Header::builder()
            .set_reference_sequences(reference_sequences)
            .build();

        print!("{}", header);
    } else {
        print!("{}", raw_header);
    }

    Ok(())
}
