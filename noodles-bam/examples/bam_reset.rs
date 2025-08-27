//! Removes alignment information from records.
//!
//! The result is similar to the output of `samtools reset <src>`. This example sets the mapping
//! quality to missing (255) rather than 0.

use std::{
    env,
    fs::File,
    io::{self, BufWriter},
};

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, io::Write, record::Flags},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = File::open(src).map(bam::io::Reader::new)?;
    let header = reader.read_header()?;

    let stdout = BufWriter::new(io::stdout().lock());
    let mut writer = sam::io::Writer::new(stdout);

    let mut new_header = header.clone();
    *new_header.header_mut() = Some(Default::default());
    new_header.reference_sequences_mut().clear();
    new_header.comments_mut().clear();

    writer.write_header(&new_header)?;

    let mut record = sam::alignment::RecordBuf::default();

    while reader.read_record_buf(&header, &mut record)? != 0 {
        let flags = record.flags();

        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        record.reference_sequence_id_mut().take();
        record.alignment_start_mut().take();
        record.mapping_quality_mut().take();
        record.cigar_mut().as_mut().clear();
        record.mate_reference_sequence_id_mut().take();
        record.mate_alignment_start_mut().take();
        *record.template_length_mut() = 0;

        let mut new_flags = Flags::UNMAPPED;

        if flags.is_segmented() {
            new_flags.insert(Flags::SEGMENTED | Flags::MATE_UNMAPPED);

            if flags.is_first_segment() {
                new_flags.insert(Flags::FIRST_SEGMENT);
            } else if flags.is_last_segment() {
                new_flags.insert(Flags::LAST_SEGMENT);
            }
        }

        *record.flags_mut() = new_flags;

        if flags.is_reverse_complemented() {
            reverse_complement(&mut record);
        }

        writer.write_alignment_record(&header, &record)?;
    }

    writer.finish(&header)?;

    Ok(())
}

fn reverse_complement(record: &mut RecordBuf) {
    fn complement(base: u8) -> u8 {
        match base {
            b'=' => b'=',
            b'A' => b'T',
            b'C' => b'G',
            b'M' => b'K',
            b'G' => b'C',
            b'R' => b'Y',
            b'S' => b'S',
            b'V' => b'B',
            b'T' => b'A',
            b'W' => b'W',
            b'Y' => b'R',
            b'H' => b'D',
            b'K' => b'M',
            b'D' => b'H',
            b'B' => b'V',
            _ => b'N',
        }
    }

    let seq = record.sequence_mut().as_mut();
    seq.reverse();

    for base in seq {
        *base = complement(*base);
    }

    let qual = record.quality_scores_mut().as_mut();
    qual.reverse();
}
