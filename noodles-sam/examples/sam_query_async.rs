//! Queries a bgzipped SAM file with a given region.
//!
//! The input bgzipped SAM file must have an associated coordinate-sorted index (CSI) in the same
//! directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, path::PathBuf};

use futures::TryStreamExt;
use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_sam as sam;
use tokio::{
    fs::File,
    io::{self, AsyncWriteExt},
};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let src = args.next().map(PathBuf::from).expect("missing src");
    let raw_region = args.next().expect("missing region");

    let mut reader = File::open(&src)
        .await
        .map(bgzf::AsyncReader::new)
        .map(sam::r#async::io::Reader::new)?;

    let header = reader.read_header().await?;

    let index = csi::r#async::read(src.with_extension("sam.csi")).await?;
    let region = raw_region.parse()?;
    let mut query = reader.query(&header, &index, &region)?;

    let mut writer = sam::r#async::io::Writer::new(io::stdout());

    while let Some(record) = query.try_next().await? {
        writer.write_alignment_record(&header, &record).await?;
    }

    writer.get_mut().shutdown().await?;

    Ok(())
}
