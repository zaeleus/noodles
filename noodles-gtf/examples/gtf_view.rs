//! Validates and prints all records in a GTF file.
//!
//! The result matches the output of `grep --invert-match "^#" <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_gtf as gtf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(gtf::io::Reader::new)?;

    let stdout = io::stdout().lock();
    let mut writer = gtf::io::Writer::new(BufWriter::new(stdout));

    let mut line = gtf::Line::default();

    while reader.read_line(&mut line)? != 0 {
        if let Some(result) = line.as_record() {
            let record = result?;
            writer.write_feature_record(&record)?;
        }
    }

    Ok(())
}
