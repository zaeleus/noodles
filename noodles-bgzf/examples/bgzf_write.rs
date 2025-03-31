//! Compresses a file as a blocked gzip file (BGZF).
//!
//! The result is similar to the output of `bgzip --stdout <src>`.

use std::{
    env,
    fs::File,
    io::{self, Write},
};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)?;

    let stdout = io::stdout().lock();
    let mut writer = bgzf::io::Writer::new(stdout);

    io::copy(&mut reader, &mut writer)?;

    writer.flush()?;

    Ok(())
}
