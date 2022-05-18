//! Downloads variants using the htsget protocol.
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

    let client = htsget::Client::new(base_url);

    let mut request = client.variants(id);

    for arg in args {
        let region = arg.parse()?;
        request = request.add_region(region);
    }

    let variants = request.send().await?;

    let mut chunks = variants.chunks();
    let mut stdout = io::stdout();

    while let Some(chunk) = chunks.try_next().await? {
        stdout.write_all(&chunk).await?;
    }

    Ok(())
}
