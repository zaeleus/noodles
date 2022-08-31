//! Prints the compressed-uncompressed position pairs of a GZ index.

use std::{env, io};

use noodles_bgzf::gzi;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let index = gzi::read(src)?;

    for record in &index {
        println!("{}\t{}", record.0, record.1);
    }

    Ok(())
}
