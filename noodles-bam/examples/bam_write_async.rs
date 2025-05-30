//! Creates a new BAM file.
//!
//! This writes a SAM header, reference sequences, and one unmapped record to stdout.
//!
//! Verify the output by piping to `samtools view --no-PG --with-header`.

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::RecordBuf,
    header::record::value::{Map, map::Program},
};
use tokio::io;

#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = bam::r#async::io::Writer::new(io::stdout());

    let header = sam::Header::builder()
        .set_header(Default::default())
        .add_program("noodles-bam", Map::<Program>::default())
        .add_comment("an example BAM written by noodles-bam")
        .build();

    writer.write_header(&header).await?;

    let record = RecordBuf::default();
    writer.write_alignment_record(&header, &record).await?;

    writer.shutdown().await?;

    Ok(())
}
