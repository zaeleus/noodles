//! Downloads a sequence using the refget protocol.
//!
//! The output is written to stdout.

use std::{
    env,
    io::{self, Write},
    ops::Bound,
};

use noodles_core::Position;
use noodles_refget as refget;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);

    let base_url = args.next().expect("missing base URL").parse()?;
    let id = args.next().expect("missing ID");

    let client = refget::Client::new(base_url);

    let mut request = client.sequence(id);

    if let Some(raw_start) = args.next() {
        let start_bound: Bound<Position> = raw_start.parse().map(Bound::Included)?;

        let end_bound = args
            .next()
            .map(|s| s.parse().map(Bound::Included))
            .transpose()?
            .unwrap_or(Bound::Unbounded);

        let interval = (start_bound, end_bound);
        request = request.set_interval(interval);
    }

    let sequence = request.send().await?;

    let stdout = io::stdout();
    let mut writer = stdout.lock();

    writer.write_all(&sequence.sequence())?;
    writeln!(writer)?;

    Ok(())
}
