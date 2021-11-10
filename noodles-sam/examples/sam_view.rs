//! Validates and prints the records in a SAM file.
//!
//! The result matches the output of `samtools view <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    reader.read_header()?;

    for result in reader.records() {
        let record = result?;
        println!("{}", record);
    }

    Ok(())
}
