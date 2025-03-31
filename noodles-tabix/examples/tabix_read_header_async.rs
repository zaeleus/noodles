//! Prints the header of the file associated with the index.
//!
//! The results match the output of `tabix --only-header <src>`.

use std::env;

use noodles_bgzf as bgzf;
use noodles_csi::BinningIndex;
use noodles_tabix as tabix;
use tokio::{
    fs::File,
    io::{self, AsyncBufReadExt},
};

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{src}.tbi");
    let index = tabix::r#async::read(tabix_src).await?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tabix header"))?;

    let reader = File::open(src).await.map(bgzf::r#async::io::Reader::new)?;
    let line_comment_prefix = char::from(header.line_comment_prefix());

    let mut lines = reader.lines();

    while let Some(line) = lines.next_line().await? {
        if !line.starts_with(line_comment_prefix) {
            break;
        }

        println!("{line}");
    }

    Ok(())
}
