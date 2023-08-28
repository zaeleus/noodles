//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --threads $(nproc) --decompress --stdout <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter},
    num::NonZeroUsize,
    thread,
};

use noodles_bgzf as bgzf;

const DEFAULT_WORKER_COUNT: NonZeroUsize = match NonZeroUsize::new(1) {
    Some(n) => n,
    None => unreachable!(),
};

fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let worker_count = args
        .next()
        .map(|s| s.parse().expect("invalid worker_count"))
        .unwrap_or_else(|| thread::available_parallelism().unwrap_or(DEFAULT_WORKER_COUNT));

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(|file| bgzf::MultithreadedReader::with_worker_count(worker_count, file))?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    io::copy(&mut reader, &mut writer)?;

    Ok(())
}
