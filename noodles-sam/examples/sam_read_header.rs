//! Prints the header of a SAM file.
//!
//! The result matches the output of `samtools head <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_sam as sam;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let mut reader = File::open(src).map(BufReader::new).map(sam::Reader::new)?;

    let header = reader.read_header()?;
    print!("{}", header);

    Ok(())
}
