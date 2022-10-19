//! Queries a bgzipped SAM file with a given region.
//!
//! The input bgzipped SAM file must have an associated coordinate-sorted index (CSI) in the same
//! directory.
//!
//! The result matches the output of `samtools view <src> <region>`.

use std::{env, fs::File, io, path::PathBuf};

use noodles_bgzf as bgzf;
use noodles_csi as csi;
use noodles_sam as sam;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args();

    let src = args.nth(1).map(PathBuf::from).expect("missing src");
    let region = args.next().expect("missing region").parse()?;

    let mut reader = File::open(&src)
        .map(bgzf::Reader::new)
        .map(sam::Reader::new)?;

    let header = reader.read_header()?.parse()?;

    let index = csi::read(src.with_extension("gz.csi"))?;
    let query = reader.query(&header, &index, &region)?;

    let stdout = io::stdout().lock();
    let mut writer = sam::Writer::new(stdout);

    for result in query {
        let record = result?;
        writer.write_record(&header, &record)?;
    }

    Ok(())
}
