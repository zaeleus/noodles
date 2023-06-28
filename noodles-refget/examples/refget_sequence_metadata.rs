//! Queries the metadata for a sequence.

use std::env;

use noodles_refget as refget;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let base_url = args.next().expect("missing base URL").parse()?;
    let id = args.next().expect("missing ID");

    let client = refget::Client::new(base_url);
    let metadata = client.sequence_metadata(id).send().await?;

    dbg!(metadata);

    Ok(())
}
