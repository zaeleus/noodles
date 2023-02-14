//! Prints the compressed-uncompressed position pairs of a GZ index.

use std::{env, io};

use noodles_bgzf::gzi;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");
    let index = gzi::read(src)?;

    for (compressed_position, uncompressed_position) in &index {
        println!("{compressed_position}\t{uncompressed_position}");
    }

    Ok(())
}
