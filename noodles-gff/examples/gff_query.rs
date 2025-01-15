//! Queries a bgzip-compressed GFF file with a given region.
//!
//! The input must have an associated coordinate-sorted index (CSI) in the same directory.
//!
//! The result matches the output of `tabix <src> <region>`.

use std::{env, fs::File, io};

use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_gff as gff;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let src = args.next().expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let index_src = format!("{src}.csi");
    let index = csi::fs::read(index_src)?;

    let mut reader = File::open(src)
        .map(bgzf::Reader::new)
        .map(gff::io::Reader::new)?;

    let query = reader.query(&index, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = gff::io::Writer::new(stdout);

    for result in query {
        let record = result?;
        writer.write_record(&record)?;
    }

    Ok(())
}
