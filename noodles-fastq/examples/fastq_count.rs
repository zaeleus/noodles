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
        .map(fastq::Reader::new)?;

    let mut record = fastq::Record::default();
    let mut n = 0;

    loop {
        match reader.read_record(&mut record) {
            Ok(0) => break,
            Ok(_) => n += 1,
            Err(e) => return Err(e),
        }
    }

    println!("{}", n);

    Ok(())
}
