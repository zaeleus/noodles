//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --decompress --stdout <src>`.

use std::{env, io};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bgzf::reader::Builder::default().build_from_path(src)?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    io::copy(&mut reader, &mut writer)?;

    Ok(())
}
