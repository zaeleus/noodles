//! Removes alignment information from records.
//!
//! The result is similar to the output of `samtools reset <src>`. This example sets the mapping
//! quality to missing (255) rather than 0.

use std::{
    env,
    io::{self, BufWriter},
};

use noodles_bam as bam;
use noodles_sam::{
    self as sam,
    alignment::RecordBuf,
    record::{sequence::Base, Flags},
    AlignmentWriter,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let src = env::args().nth(1).expect("missing src");

    let mut reader = bam::reader::Builder.build_from_path(src)?;
    let header = reader.read_header()?;

    let stdout = BufWriter::new(io::stdout().lock());
    let mut writer = sam::Writer::new(stdout);

    let mut new_header = header.clone();
    *new_header.header_mut() = Some(Default::default());
    new_header.reference_sequences_mut().clear();
    new_header.comments_mut().clear();

    writer.write_header(&new_header)?;

    let mut record = sam::alignment::RecordBuf::default();

    while reader.read_record(&header, &mut record)? != 0 {
        let flags = record.flags();

        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        record.reference_sequence_id_mut().take();
        record.alignment_start_mut().take();
        record.mapping_quality_mut().take();
        record.cigar_mut().clear();
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
    fn complement(base: Base) -> Base {
        match base {
            Base::Eq => Base::Eq,
            Base::A => Base::T,
            Base::C => Base::G,
            Base::M => Base::K,
            Base::G => Base::C,
            Base::R => Base::Y,
            Base::S => Base::S,
            Base::V => Base::B,
            Base::T => Base::A,
            Base::W => Base::W,
            Base::Y => Base::R,
            Base::H => Base::D,
            Base::K => Base::M,
            Base::D => Base::H,
            Base::B => Base::V,
            _ => Base::N,
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
