//! Converts a BAM file to the FASTQ format.
//!
//! This example only supports single segment reads.
//!
//! The result matches the output of `samtools fastq <src>`.

use std::{
    env,
    fs::File,
    io::{self, BufWriter, Write},
};

use noodles_bam as bam;
use noodles_sam::alignment::record::Flags;

const FILTERS: Flags = Flags::SECONDARY.union(Flags::SUPPLEMENTARY);

fn main() -> io::Result<()> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    reader.read_header()?;

    let stdout = io::stdout().lock();
    let mut writer = BufWriter::new(stdout);

    for result in reader.records() {
        let record = result?;

        let flags = record.flags();

        if flags.is_segmented() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "this example only supports single segment reads",
            ));
        } else if flags.intersects(FILTERS) {
            continue;
        }

        write_record(&mut writer, &record)?;
    }

    Ok(())
}

fn write_record<W>(writer: &mut W, record: &bam::Record) -> io::Result<()>
where
    W: Write,
{
    const DEFINITION_PREFIX: u8 = b'@';
    const MISSING_NAME: &[u8] = b"*";
    const SEPARATOR: u8 = b'+';
    const LINE_FEED: u8 = b'\n';

    writer.write_all(&[DEFINITION_PREFIX])?;

    let name = record
        .name()
        .map(|name| name.as_ref())
        .unwrap_or(MISSING_NAME);

    writer.write_all(name)?;
    writer.write_all(&[LINE_FEED])?;

    let is_reversed_complemented = record.flags().is_reverse_complemented();

    let bases = record.sequence().iter();

    if is_reversed_complemented {
        for base in bases.rev().map(complement_base) {
            writer.write_all(&[base])?;
        }
    } else {
        for base in bases {
            writer.write_all(&[base])?;
        }
    }

    writer.write_all(&[LINE_FEED])?;

    writer.write_all(&[SEPARATOR, LINE_FEED])?;

    let quality_scores = record.quality_scores();
    let scores = quality_scores.as_ref().iter().copied().map(encode_score);

    if is_reversed_complemented {
        for n in scores.rev() {
            writer.write_all(&[n])?;
        }
    } else {
        for n in scores {
            writer.write_all(&[n])?;
        }
    }

    writer.write_all(&[LINE_FEED])?;

    Ok(())
}

fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'U' => b'A',
        b'W' => b'W',
        b'S' => b'S',
        b'M' => b'K',
        b'K' => b'M',
        b'R' => b'Y',
        b'Y' => b'R',
        b'B' => b'V',
        b'D' => b'H',
        b'H' => b'D',
        b'V' => b'B',
        b'N' => b'N',
        _ => unreachable!(),
    }
}

fn encode_score(n: u8) -> u8 {
    const OFFSET: u8 = b'!';
    n.checked_add(OFFSET).expect("attempt to add with overflow")
}
