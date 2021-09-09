//! Prints the header of a SAM file.
//!
//! The result matches the output of `samtools view --no-PG --header-only <src>`.

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
    print!("{}", header);

    Ok(())
}
