//! Downloads reads using the htsget protocol.
//!
//! The output is written to stdout.

use std::env;

use futures::TryStreamExt;
use noodles_htsget as htsget;
use tokio::io::{self, AsyncWriteExt};

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let base_url = args.next().expect("missing base URL").parse()?;
    let id = args.next().expect("missing ID");
    let region = args.next().map(|s| s.parse()).transpose()?;

    let client = htsget::Client::new(base_url);

    let mut request = client.reads(id);

    if let Some(region) = region {
        request = request.set_region(region);
    }

    let reads = request.send().await?;

    let mut chunks = reads.chunks();
    let mut stdout = io::stdout();

    while let Some(chunk) = chunks.try_next().await? {
        stdout.write_all(&chunk).await?;
    }

    Ok(())
}
