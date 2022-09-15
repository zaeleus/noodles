//! Seeks and decompresses a BGZF file starting at an uncompressed position and ending in `n`
//! bytes.
//!
//! The source file must have an associated index.
//!
//! The result matches the output of `bgzip --offset <position> --size <length> <src>`.

use std::{
    env,
    fs::File,
    io::{self, Read, Write},
};

use noodles_bgzf::{self as bgzf, gzi};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let raw_position = args.next().expect("missing position");
    let raw_length = args.next().expect("missing length");

    let pos = raw_position.parse()?;
    let len = raw_length.parse()?;

    let mut reader = File::open(&src).map(bgzf::Reader::new)?;
    let index = gzi::read(format!("{}.gzi", src))?;

    reader.seek_by_uncompressed_position(&index, pos)?;

    let mut buf = vec![0; len];
    reader.read_exact(&mut buf)?;

    let stdout = io::stdout();
    let mut handle = stdout.lock();
    handle.write_all(&buf)?;

    Ok(())
}
