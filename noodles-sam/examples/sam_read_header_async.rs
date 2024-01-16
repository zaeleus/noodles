//! Prints the header of a SAM file.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::env;

use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, BufReader},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(sam::AsyncReader::new)?;

    let header = reader.read_header().await?;

    let mut writer = sam::AsyncWriter::new(io::stdout());
    writer.write_header(&header).await?;

    Ok(())
}
