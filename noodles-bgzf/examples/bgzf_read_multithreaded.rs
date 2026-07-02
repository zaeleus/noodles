//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --threads $(nproc) --decompress --stdout <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(bgzf::io::MultithreadedReader::new)?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    io::copy(&mut reader, &mut writer)?;

    Ok(())
}
