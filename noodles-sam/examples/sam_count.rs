//! Counts the number of records in a SAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::{env, fs::File, io::BufReader};

use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    let header = reader.read_header()?.parse()?;

    let mut n = 0;

    for result in reader.records(&header) {
        let _ = result?;
        n += 1;
    }

    println!("{}", n);

    Ok(())
}
