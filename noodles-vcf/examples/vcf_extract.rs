//! Extracts fields from a VCF.
//!
//! The output is similar to `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n" <src>`.

use std::{
    env,
    io::{self, BufWriter, Write},
};

use noodles_vcf::{
    self as vcf,
    variant::record_buf::samples::{keys::key, sample::Value},
};

const MISSING: &str = ".";

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    for result in reader.record_bufs(&header) {
        let record = result?;

        write!(
            writer,
            "{chrom}\t{pos}\t{ref}",
            chrom = record.chromosome(),
            pos = record.position(),
            ref = record.reference_bases(),
        )?;

        write!(writer, "\t")?;

        if record.alternate_bases().is_empty() {
            write!(writer, "{MISSING}")?;
        } else {
            for (i, allele) in record.alternate_bases().as_ref().iter().enumerate() {
                if i > 0 {
                    write!(writer, ",")?;
                }

                write!(writer, "{allele}")?;
            }
        }

        for (sample_name, sample) in header.sample_names().iter().zip(record.samples().values()) {
            let value = sample
                .get(key::GENOTYPE)
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
