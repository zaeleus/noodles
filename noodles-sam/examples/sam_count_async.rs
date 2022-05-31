//! Counts the number of records in a SAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::env;

use futures::TryStreamExt;
use noodles_sam as sam;
use tokio::{fs::File, io::BufReader};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .await
        .map(BufReader::new)
        .map(sam::AsyncReader::new)?;

    let header = reader.read_header().await?.parse()?;

    let mut records = reader.records(&header);
    let mut n = 0;

    while records.try_next().await?.is_some() {
        n += 1;
    }

    println!("{}", n);

    Ok(())
}
