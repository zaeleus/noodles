//! Validates and prints the records in a SAM file.
//!
//! The result matches the output of `samtools view <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = sam::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(BufWriter::new(stdout));

    for result in reader.record_bufs(&header) {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
