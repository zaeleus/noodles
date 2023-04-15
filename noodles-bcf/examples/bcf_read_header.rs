//! Prints the header of a BCF file.
//!
//! The result matches the output of `bcftools head <src>`.

use std::{env, fs::File, io};

use noodles_bcf as bcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bcf::Reader::new)?;
    reader.read_file_format()?;

    let header = reader.read_header()?;
    print!("{header}");

    Ok(())
}
