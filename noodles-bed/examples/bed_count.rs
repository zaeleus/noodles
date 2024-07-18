//! Counts the number of records in a BED3+ file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_bed as bed;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(bed::io::Reader::<_, 3>::new)?;

    let mut record = bed::Record::default();
    let mut n = 0;

    while reader.read_record(&mut record)? != 0 {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
