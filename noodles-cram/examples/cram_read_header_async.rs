//! Prints the header of a CRAM file.
//!
//! A CRAM file header is a SAM header.
//!
//! The result is similar to the output of `samtools head <src>`.

use std::env;

use noodles_cram as cram;
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(cram::AsyncReader::new)?;
    reader.read_file_definition().await?;

    let header = reader.read_file_header().await?;
    print!("{header}");

    Ok(())
}
