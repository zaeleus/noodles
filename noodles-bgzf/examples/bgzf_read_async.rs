//! Decompresses a BGZF file using an asynchronous reader.
//!
//! The result matches the output of `bgzip --decompress --stdout <src>`.

use std::env;

use noodles_bgzf as bgzf;
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bgzf::r#async::fs::open(src).await?;
    let mut writer = io::BufWriter::new(io::stdout());
    io::copy_buf(&mut reader, &mut writer).await?;

    Ok(())
}
