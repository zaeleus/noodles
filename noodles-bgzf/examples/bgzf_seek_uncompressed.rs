//! Seeks and decompresses a BGZF file starting at an uncompressed position and ending in `n`
//! bytes.
//!
//! The source file must have an associated index.
//!
//! The result matches the output of `bgzip --offset <position> --size <length> <src>`.

use std::{
    env,
    io::{self, Read, Seek, SeekFrom, Write},
};

use noodles_bgzf as bgzf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let position = args.next().expect("missing position").parse()?;
    let length = args.next().expect("missing length").parse()?;

    let mut reader = bgzf::io::indexed_reader::Builder::default().build_from_path(src)?;
    reader.seek(SeekFrom::Start(position))?;

    let mut buf = vec![0; length];
    reader.read_exact(&mut buf)?;

    let mut writer = io::stdout().lock();
    writer.write_all(&buf)?;

    Ok(())
}
