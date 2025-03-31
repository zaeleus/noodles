//! Compresses a file in the blocked gzip format (BGZF).
//!
//! The result is similar to the output of `bgzip --stdout <src>`.

use std::env;

use noodles_bgzf as bgzf;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await?;
    let mut writer = bgzf::r#async::io::Writer::new(io::stdout());

    io::copy(&mut reader, &mut writer).await?;

    writer.shutdown().await?;

    Ok(())
}
