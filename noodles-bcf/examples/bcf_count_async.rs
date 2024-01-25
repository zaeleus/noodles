//! Counts the number of records in a BCF file.
//!
//! The result matches the output of `bcftools view --no-header <src> | wc -l`.

use std::env;

use futures::TryStreamExt;
use noodles_bcf as bcf;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(bcf::r#async::io::Reader::new)?;
    reader.read_file_format().await?;
    reader.read_header().await?;

    let mut records = reader.lazy_records();
    let mut n = 0;

    while records.try_next().await?.is_some() {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
