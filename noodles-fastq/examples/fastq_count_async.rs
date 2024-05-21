//! Counts the number of records in a FASTQ file.

use std::env;

use futures::TryStreamExt;
use noodles_fastq as fastq;
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
        .map(fastq::r#async::io::Reader::new)?;

    let mut records = reader.records();
    let mut n = 0;

    while records.try_next().await?.is_some() {
        n += 1;
    }

    println!("{n}");

    Ok(())
}
