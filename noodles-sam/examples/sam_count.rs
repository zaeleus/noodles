//! Counts the number of records in a SAM file.
//!
//! The result matches the output of `samtools view -c <src>`.

use std::{env, fs::File, io::BufReader};

use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;
    reader.read_header()?;

    let mut buf = String::new();
    let mut n = 0;

    loop {
        buf.clear();

        let bytes_read = reader.read_record(&mut buf)?;

        if bytes_read == 0 {
            break;
        }

        n += 1;
    }

    println!("{}", n);

    Ok(())
}
