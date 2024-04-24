//! Compresses a file as a blocked gzip file (BGZF).
//!
//! The result is similar to the output of `bgzip --threads $(nproc) --stdout <src>`.

use std::{env, fs::File, io, num::NonZeroUsize, thread};

use noodles_bgzf as bgzf;

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let worker_count = args
        .next()
        .map(|s| s.parse().expect("invalid worker_count"))
        .unwrap_or_else(|| thread::available_parallelism().unwrap_or(NonZeroUsize::MIN));

    let mut reader = File::open(src)?;

    let mut writer = bgzf::MultithreadedWriter::with_worker_count(worker_count, io::stdout());
    io::copy(&mut reader, &mut writer)?;
    writer.finish()?;

    Ok(())
}
