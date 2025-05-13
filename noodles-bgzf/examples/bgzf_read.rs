//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --decompress --stdout <src>`.

use std::{env, io};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bgzf::fs::open(src)?;
    let mut writer = io::stdout().lock();
    io::copy(&mut reader, &mut writer)?;

    Ok(())
}
