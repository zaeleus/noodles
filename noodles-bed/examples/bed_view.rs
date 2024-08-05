//! Prints BED records as BED3+ records.

use std::{env, io};

use noodles_bed as bed;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bed::io::reader::Builder::<3>.build_from_path(src)?;

    let stdout = io::stdout().lock();
    let mut writer = bed::io::Writer::<_, 3>::new(stdout);

    let mut record = bed::Record::default();

    while reader.read_record(&mut record)? != 0 {
        writer.write_record(&record)?;
    }

    Ok(())
}
