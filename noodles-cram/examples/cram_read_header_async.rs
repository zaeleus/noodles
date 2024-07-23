//! Prints the header of a CRAM file.
//!
//! A CRAM file header is a SAM header.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::env;

use noodles_cram as cram;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(cram::AsyncReader::new)?;
    let header = reader.read_header().await?;

    let mut writer = sam::AsyncWriter::new(io::stdout());
    writer.write_header(&header).await?;

    writer.get_mut().shutdown().await?;

    Ok(())
}
