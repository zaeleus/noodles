use std::fs::File;
use std::io;
use std::path::Path;

use formats::{bam, fastq};

pub fn convert<P, Q>(src: P, dst: Q) -> io::Result<()>
where
    P: AsRef<Path>,
    Q: AsRef<Path>,
{
    let mut reader = bam::Reader::<File>::open(src)?;
    let mut writer = fastq::Writer::<File>::create(dst)?;

    let _header = reader.header()?;
    let _references = reader.references()?.count();

    for record in reader.records() {
        writer.write_bam_record(&record)?;
    }

    Ok(())
}
