//! Counts the number of records in a SAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(sam::io::Reader::new)?;

    reader.read_header()?;

    let mut n = 0;

    for result in reader.records() {
        let _ = result?;
        n += 1;
    }

    println!("{n}");

    Ok(())
}
