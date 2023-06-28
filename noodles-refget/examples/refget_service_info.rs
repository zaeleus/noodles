//! Returns refget service information.

use std::env;

use noodles_refget as refget;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base_url = env::args().nth(1).expect("missing base URL").parse()?;

    let client = refget::Client::new(base_url);
    let service = client.service_info().send().await?;

    dbg!(service);

    Ok(())
}
