//! Counts the number of records in a BED3+ file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_bed as bed;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(bed::Reader::new)?;

    let mut buf = String::new();
    let mut n = 0;

    while reader.read_record(&mut buf)? != 0 {
        buf.parse::<bed::Record<3>>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        buf.clear();
        n += 1;
    }

    println!("{}", n);

    Ok(())
}
