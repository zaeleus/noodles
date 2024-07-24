//! Rewrites a variant format to another variant format.
//!
//! The output format is determined from the extension of the destination.

use std::env;

use futures::TryStreamExt;
use noodles_util::variant;
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let mut args = env::args().skip(1);

    let src = args.next().expect("missing src");
    let dst = args.next().expect("missing dst");

    let mut reader = variant::r#async::io::reader::Builder::default()
        .build_from_path(src)
        .await?;

    let header = reader.read_header().await?;

    let mut writer = variant::r#async::io::writer::Builder::default()
        .build_from_path(dst)
        .await?;

    writer.write_header(&header).await?;

    let mut records = reader.records();

    while let Some(record) = records.try_next().await? {
        writer.write_record(&header, record.as_ref()).await?;
    }

    Ok(())
}
