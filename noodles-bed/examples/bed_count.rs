//! Counts the number of records in a BED3+ file.

use std::{env, io};

use noodles_bed as bed;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bed::io::reader::Builder::<3>.build_from_path(src)?;

    let mut record = bed::Record::default();
    let mut n = 0;

    while reader.read_record(&mut record)? != 0 {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
