//! Prints the reference sequence names stored in the index.
//!
//! The results match the output of `tabix --list-chroms <src>`.

use std::env;

use noodles_csi::BinningIndex;
use noodles_tabix as tabix;
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{src}.tbi");
    let index = tabix::r#async::read(tabix_src).await?;

    let header = index
        .header()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing tabix header"))?;

    for reference_sequence_name in header.reference_sequence_names() {
        println!("{reference_sequence_name}");
    }

    Ok(())
}
