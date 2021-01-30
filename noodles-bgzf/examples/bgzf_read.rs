//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --decompress --stdout <src>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bgzf::Reader::new)?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    io::copy(&mut reader, &mut writer)?;

    Ok(())
}
