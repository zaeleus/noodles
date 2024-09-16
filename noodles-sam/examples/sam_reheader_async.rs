//! Replaces the header of a SAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::env;

use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt, BufReader},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(sam::r#async::io::Reader::new)?;

    let mut header: sam::Header = reader.read_header().await?;
    header.add_comment("a comment added by noodles-sam");

    let mut writer = sam::r#async::io::Writer::new(io::stdout());
    writer.write_header(&header).await?;

    io::copy(reader.get_mut(), writer.get_mut()).await?;

    writer.get_mut().shutdown().await?;

    Ok(())
}
