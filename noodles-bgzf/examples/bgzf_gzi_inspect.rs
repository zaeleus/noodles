//! Prints the compressed-uncompressed position pairs of a GZ index.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bgzf::gzi::Reader::new)?;
    let index = reader.read_index()?;

    for record in &index {
        println!("{}\t{}", record.0, record.1);
    }

    Ok(())
}
