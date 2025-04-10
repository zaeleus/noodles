//! Creates a new GTF file.
//!
//! This writes one (sparse) record to stdout.

use std::io;

use noodles_gff as gff;
use noodles_gtf as gtf;

fn main() -> io::Result<()> {
    let stdout = io::stdout().lock();
    let mut writer = gtf::io::Writer::new(stdout);

    let record = gff::feature::RecordBuf::default();
    writer.write_record(&record)?;

    Ok(())
}
