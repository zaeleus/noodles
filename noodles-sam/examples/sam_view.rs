//! Validates and prints the records in a SAM file.
//!
//! The result matches the output of `samtools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    let header = reader.read_header()?.parse()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(BufWriter::new(stdout));

    for result in reader.records(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
