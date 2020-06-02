//! Prints the header of a BAM file.
//!
//! A BAM file header is a SAM header.
//!
//! The result matches the output of `samtools view -H --no-PG <src>`.

use std::{env, fs::File, io};

use noodles_bam as bam;
use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::Reader::new)?;
    let header = reader.read_header()?;

    if header.is_empty() {
        let reference_sequences = reader.read_reference_sequences()?;
        let mut builder = sam::Header::builder();

        for reference_sequence in reference_sequences {
            builder = builder.add_reference_sequence(reference_sequence);
        }

        print!("{}", builder.build());
    } else {
        print!("{}", header);
    }

    Ok(())
}
