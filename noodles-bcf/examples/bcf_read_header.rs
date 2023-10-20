//! Prints the header of a BCF file.
//!
//! The result matches the output of `bcftools head <src>`.

use std::{env, io};

use noodles_bcf as bcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bcf::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;
    print!("{header}");

    Ok(())
}
