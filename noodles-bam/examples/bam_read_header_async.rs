//! Prints the header of a BAM file.
//!
//! A BAM file header is a SAM header.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::env;

use noodles_bam as bam;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bam::r#async::io::Reader::new)?;
    let header = reader.read_header().await?;

    let mut writer = sam::r#async::io::Writer::new(io::stdout());
    writer.write_header(&header).await?;

    writer.get_mut().shutdown().await?;

    Ok(())
}
