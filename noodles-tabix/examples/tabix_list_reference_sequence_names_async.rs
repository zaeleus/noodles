//! Prints the reference sequence names stored in the index.
//!
//! The results match the output of `tabix --list-chroms <src>`.

use std::env;

use noodles_tabix as tabix;
use tokio::io;

#[tokio::main]
async fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let tabix_src = format!("{}.tbi", src);
    let index = tabix::r#async::read(tabix_src).await?;

    for reference_sequence_name in index.header().reference_sequence_names() {
        println!("{}", reference_sequence_name);
    }

    Ok(())
}
