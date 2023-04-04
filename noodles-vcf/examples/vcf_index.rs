//! Builds and write a tabix index from a bgzipped VCF file.
//!
//! This writes the output to stdout rather than `<src>.tbi`.
//!
//! The output is similar to the output of `bcftools index --tbi <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::{self as csi, index::reference_sequence::bin::Chunk};
use noodles_tabix as tabix;
use noodles_vcf as vcf;

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(bgzf::Reader::new)
        .map(vcf::Reader::new)?;

    let header = reader.read_header()?;

    let mut record = vcf::Record::default();

    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::index::header::Builder::vcf().build());

    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&header, &mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let reference_sequence_name = record.chromosome().to_string();
        let start = Position::try_from(usize::from(record.position()))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let end = record
            .end()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|position| {
                Position::try_from(usize::from(position))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            })?;

        indexer.add_record(&reference_sequence_name, start, end, chunk)?;

        start_position = end_position;
    }

    let index = indexer.build();

    let stdout = io::stdout().lock();
    let mut writer = tabix::Writer::new(BufWriter::new(stdout));

    writer.write_index(&index)?;

    Ok(())
}
