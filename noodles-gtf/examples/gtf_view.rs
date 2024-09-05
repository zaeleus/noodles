//! Validates and prints all records in a GTF file.
//!
//! The result matches the output of `grep --invert-match "^#" <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_gtf as gtf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(gtf::io::Reader::new)?;

    for result in reader.records() {
        let record = result?;
        println!("{record}");
    }

    Ok(())
}
