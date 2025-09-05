//! Extracts fields from a VCF.
//!
//! The output is similar to `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n" <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufReader, BufWriter, Write},
};

use noodles_vcf::{
    self as vcf,
    variant::record::{
        AlternateBases,
        samples::{Sample, keys::key, series::Value},
    },
};

const MISSING: &str = ".";

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;

    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    for result in reader.records() {
        let record = result?;

        let position = record.variant_start().transpose()?;

        write!(
            writer,
            "{chrom}\t{pos}\t{ref}",
            chrom = record.reference_sequence_name(),
            pos = position.map(usize::from).unwrap_or_default(),
            ref = record.reference_bases(),
        )?;

        write!(writer, "\t")?;

        if record.alternate_bases().is_empty() {
            write!(writer, "{MISSING}")?;
        } else {
            for (i, result) in record.alternate_bases().iter().enumerate() {
                if i > 0 {
                    write!(writer, ",")?;
                }

                let allele = result?;
                write!(writer, "{allele}")?;
            }
        }

        for (sample_name, sample) in header.sample_names().iter().zip(record.samples().iter()) {
            let value = sample
                .get(&header, key::GENOTYPE)
                .transpose()?
                .expect("missing GT field")
                .expect("missing GT value");

            write!(writer, "\t{sample_name}=")?;

            match value {
                Value::String(gt) => write!(writer, "{gt}")?,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("expected GT to be a string, got {value:?}"),
                    ));
                }
            }
        }

        writeln!(writer)?;
    }

    Ok(())
}
