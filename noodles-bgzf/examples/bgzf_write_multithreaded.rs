//! Compresses a file as a blocked gzip file (BGZF).
//!
//! The result is similar to the output of `bgzip --threads $(nproc) --stdout <src>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");

    let mut reader = File::open(src)?;

    let mut writer = bgzf::io::MultithreadedWriter::new(io::stdout());
    io::copy(&mut reader, &mut writer)?;
    writer.finish()?;

    Ok(())
}
