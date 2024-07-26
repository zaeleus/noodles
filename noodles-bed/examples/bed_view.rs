//! Prints BED records as BED3+ records.

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
        .map(bed::io::Reader::<_, 3>::new)?;

    let stdout = io::stdout().lock();
    let mut writer = bed::io::Writer::new(stdout);

    let mut record = bed::Record::default();

    while reader.read_record(&mut record)? != 0 {
        writer.write_record(&record)?;
    }

    Ok(())
}
