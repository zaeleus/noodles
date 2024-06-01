//! Validates and prints the records as BED3+ records in a BED file.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_bed as bed;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(bed::io::Reader::new)?;

    let stdout = io::stdout().lock();
    let mut writer = bed::io::Writer::new(stdout);

    for result in reader.records::<3>() {
        let record = result?;
        writer.write_record(&record)?;
    }

    Ok(())
}
