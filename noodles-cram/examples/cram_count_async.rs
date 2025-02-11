//! Counts the number of records in a CRAM file.
//!
//! Note that this example counts each successfully read record and not by summing the number of
//! records field in each container or slice header.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::env;

use noodles_cram as cram;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(cram::r#async::io::Reader::new)?;
    reader.read_header().await?;

    let mut n = 0;

    while let Some(container) = reader.read_container().await? {
        for slice in container.slices() {
            let records = slice.records(container.compression_header())?;
            n += records.len();
        }
    }

    println!("{n}");

    Ok(())
}
