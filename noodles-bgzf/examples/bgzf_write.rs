//! Compresses a file as a blocked gzip file (BGZF).
//!
//! The result is similar to the output of `bgzip --stdout <src>`.

use std::{
    env,
    fs::File,
    io::{self, Read, Write},
};

use noodles_bgzf as bgzf;

const DEFAULT_BUF_SIZE: usize = 8192;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)?;

    let stdout = io::stdout();
    let handle = stdout.lock();
    let mut writer = bgzf::Writer::new(handle);

    let mut buf = vec![0; DEFAULT_BUF_SIZE];

    loop {
        let bytes_read = match reader.read(&mut buf) {
            Ok(0) => break,
            Ok(n) => n,
            Err(e) => return Err(e),
        };

        writer.write_all(&buf[..bytes_read])?;
    }

    writer.finish()?;

    Ok(())
}
