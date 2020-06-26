//! Counts the number of records in a GFF file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_gff as gff;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(gff::Reader::new)?;

    let count = reader.records().count();
    println!("{}", count);

    Ok(())
}
