//! Validates and prints the records as BED3+ records in a BED file.

use std::{env, fs::File, io::BufReader};

use noodles_bed as bed;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(bed::io::Reader::new)?;

    for result in reader.records::<3>() {
        let record = result?;
        println!("{record}");
    }

    Ok(())
}
