//! Counts the number of records in a BCF file.
//!
//! The result matches the output of `bcftools view --no-header <src> | wc -l`.

use std::{env, fs::File, io};

use noodles_bcf as bcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bcf::Reader::new)?;
    reader.read_file_format()?;
    reader.read_header()?;

    let mut n = 0;

    for result in reader.lazy_records() {
        let _ = result?;
        n += 1;
    }

    println!("{n}");

    Ok(())
}
