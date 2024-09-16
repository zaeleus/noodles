//! Replaces the SAM header of a BAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::env;

use noodles_bam as bam;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::r#async::io::Reader::new)?;
    let mut header = reader.read_header().await?;

    header.add_comment("a comment added by noodles-bam");

    let mut writer = bam::r#async::io::Writer::new(io::stdout());
    writer.write_header(&header).await?;

    io::copy(reader.get_mut(), writer.get_mut()).await?;

    writer.shutdown().await?;

    Ok(())
}
