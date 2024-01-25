//! Extracts fields from a VCF.
//!
//! The output is similar to `bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n" <src>`.

use std::{
    env,
    io::{self, BufWriter, Write},
};

use noodles_vcf::{self as vcf, record::genotypes::keys::key};

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = vcf::io::reader::Builder::default().build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    for result in reader.records(&header) {
        let record = result?;

        write!(
            writer,
            "{chrom}\t{pos}\t{ref}\t{alt}",
            chrom = record.chromosome(),
            pos = record.position(),
            ref = record.reference_bases(),
            alt = record.alternate_bases(),
        )?;

        for (sample_name, sample) in header
            .sample_names()
            .iter()
            .zip(record.genotypes().values())
        {
            let gt = sample
                .get(&key::GENOTYPE)
                .expect("missing GT field")
                .expect("missing GT value");

            write!(writer, "\t{sample_name}={gt}")?;
        }

        writeln!(writer)?;
    }

    Ok(())
}
