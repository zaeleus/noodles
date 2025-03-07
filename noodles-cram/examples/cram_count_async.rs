//! Counts the number of records in a CRAM file.
//!
//! The result matches the output of `samtools view --count <src>`.

use std::env;

use noodles_cram::{self as cram, io::reader::Container};
use tokio::{fs::File, io};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).await.map(cram::r#async::io::Reader::new)?;
    reader.read_header().await?;

    let mut container = Container::default();
    let mut n = 0;

    while reader.read_container(&mut container).await? != 0 {
        n += container.header().record_count();
    }

    println!("{n}");

    Ok(())
}
