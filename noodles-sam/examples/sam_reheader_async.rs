//! Replaces the header of a SAM file.
//!
//! This is similar to the functionality of `samtools reheader`.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use std::env;

use futures::TryStreamExt;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, BufReader},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(sam::AsyncReader::new)?;

    let mut header: sam::Header = reader.read_header().await?.parse()?;
    header.add_comment("a comment added by noodles-sam");

    let mut writer = sam::AsyncWriter::new(io::stdout());
    writer.write_header(&header).await?;

    let mut records = reader.records(&header);

    while let Some(record) = records.try_next().await? {
        writer.write_record(&header, &record).await?;
    }

    Ok(())
}
