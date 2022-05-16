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
    let reference_name = args.next().expect("missing reference name");
    let start = args.next().expect("missing start").parse()?;
    let end = args.next().expect("missing end").parse()?;

    let client = htsget::Client::new(base_url);

    let reads = client
        .variants(id)
        .set_reference_name(reference_name)
        .set_start(start)
        .set_end(end)
        .send()
        .await?;

    let mut chunks = reads.chunks();
    let mut stdout = io::stdout();

    while let Some(chunk) = chunks.try_next().await? {
        stdout.write_all(&chunk).await?;
    }

    Ok(())
}
