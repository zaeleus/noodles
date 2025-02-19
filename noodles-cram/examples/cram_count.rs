//! Counts the number of records in a CRAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::{env, io};

use noodles_cram::{self as cram, io::reader::Container};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = cram::io::reader::Builder::default().build_from_path(src)?;
    reader.read_header()?;

    let mut container = Container::default();
    let mut n = 0;

    while reader.read_container(&mut container)? != 0 {
        n += container.header().record_count();
    }

    println!("{n}");

    Ok(())
}
