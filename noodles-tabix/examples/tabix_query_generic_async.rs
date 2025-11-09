//! Queries an indexed TSV with the given region.
//!
//! The input must have an associated tabix index in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::env;

use futures::TryStreamExt;
use noodles_csi as csi;
use noodles_tabix as tabix;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let index_src = format!("{src}.tbi");
    let index = tabix::fs::read(index_src)?;

    let mut reader = File::open(src)
        .await
        .map(|file| csi::r#async::io::IndexedReader::new(file, index))?;

    let mut query = reader.query(&region)?;

    while let Some(record) = query.try_next().await? {
        println!("{}", record.as_ref());
    }

    Ok(())
}
