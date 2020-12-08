//! Decompresses a BGZF file.
//!
//! The result matches the output of `bgzip --decompress --stdout <src>`.

use std::{
    env,
    fs::File,
    io::{self, Read, Write},
};

use noodles_bgzf as bgzf;

const DEFAULT_BUF_SIZE: usize = 8192;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bgzf::Reader::new)?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    let mut buf = vec![0; DEFAULT_BUF_SIZE];

    loop {
        match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => writer.write_all(&buf[..n])?,
            Err(ref e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }

    Ok(())
}
