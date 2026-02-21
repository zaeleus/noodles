//! Filters VCF records by quality score.
//!
//! This example demonstrates how to filter VCF records based on a minimum quality score
//! threshold. Records with a quality score below the threshold are excluded from the output.
//!
//! The result is similar to using `bcftools view -i 'QUAL>=threshold' <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader},
};

use noodles_vcf as vcf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    
    let src = args.next().expect("missing src");
    let min_quality: f32 = args
        .next()
        .expect("missing minimum quality threshold")
        .parse()
        .expect("invalid quality threshold");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = vcf::io::Writer::new(stdout);

    writer.write_header(&header)?;

    for result in reader.records() {
        let record = result?;
        
        // Filter by quality score
        if let Some(quality_score) = record.quality_score() {
            if let Some(score) = quality_score.value() {
                if score >= min_quality {
                    writer.write_record(&header, &record)?;
                }
            }
        }
    }

    Ok(())
}
