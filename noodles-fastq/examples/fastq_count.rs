//! Counts the number of records in a FASTQ file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_fastq as fastq;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(fastq::io::Reader::new)?;

    let mut n = 0;

    for result in reader.records() {
        let _ = result?;
        n += 1;
    }

    println!("{n}");

    Ok(())
}
