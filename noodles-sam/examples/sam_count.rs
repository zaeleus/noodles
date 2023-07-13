//! Counts the number of records in a SAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::{env, io};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = sam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let mut n = 0;

    for result in reader.records(&header) {
        let _ = result?;
        n += 1;
    }

    println!("{n}");

    Ok(())
}
