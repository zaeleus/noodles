//! Builds and writes a tabix index from a bgzip-compressed BED file.
//!
//! This writes the output to stdout rather than `<src>.tbi`.
//!
//! The output is similar to the output of `tabix --preset bed <src>`.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_bed as bed;
use noodles_tabix as tabix;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let index = bed::fs::index(src)?;

    let stdout = io::stdout().lock();
    let mut writer = tabix::io::Writer::new(BufWriter::new(stdout));
    writer.write_index(&index)?;

    Ok(())
}
