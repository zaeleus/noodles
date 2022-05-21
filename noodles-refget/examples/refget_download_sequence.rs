//! Downloads a sequence using the refget protocol.
//!
//! The output is written to stdout.

use std::{
    env,
    io::{self, Write},
};

use noodles_core::region::Interval;
use noodles_refget as refget;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let base_url = args.next().expect("missing base URL").parse()?;
    let id = args.next().expect("missing ID");
    let interval: Option<Interval> = args.next().map(|s| s.parse()).transpose()?;

    let client = refget::Client::new(base_url);

    let mut request = client.sequence(id);

    if let Some(interval) = interval {
        request = request.set_interval(interval);
    }

    let sequence = request.send().await?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    writer.write_all(&sequence.sequence())?;
    writeln!(writer)?;

    Ok(())
}
