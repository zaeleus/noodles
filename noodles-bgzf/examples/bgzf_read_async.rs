//! Decompresses a BGZF file using an asynchronous reader.
//!
//! The result matches the ouptput of `bgzip --decompress --stdout <src>`.

use std::env;

use noodles_bgzf as bgzf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bgzf::AsyncReader::new)?;
    let mut writer = io::BufWriter::new(io::stdout());
    io::copy_buf(&mut reader, &mut writer).await?;

    Ok(())
}
