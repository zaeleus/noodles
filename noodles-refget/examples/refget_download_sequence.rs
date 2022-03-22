//! Downloads a sequence using the refget protocol.
//!
//! The output is written to stdout.

use std::{
    env,
    io::{self, Write},
};

use noodles_refget as refget;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let base_url = args.next().expect("missing base URL").parse()?;
    let id = args.next().expect("missing ID");
    let start = args.next().expect("missing start").parse()?;
    let end = args.next().expect("missing end").parse()?;

    let client = refget::Client::new(base_url);

    let sequence = client
        .sequence(id)
        .set_start(start)
        .set_end(end)
        .send()
        .await?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    writer.write_all(&sequence.sequence())?;
    writeln!(writer)?;

    Ok(())
}
