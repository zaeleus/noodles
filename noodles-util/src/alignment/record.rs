use std::io;

use bstr::BStr;
use noodles_bam as bam;
use noodles_core::Position;
use noodles_sam as sam;

/// An alignment record.
pub enum Record {
    /// A SAM record.
    Sam(sam::Record),
    /// A BAM record.
    Bam(bam::Record),
    /// A CRAM record.
    Cram(sam::alignment::RecordBuf),
}

impl sam::alignment::Record for Record {
    fn name(&self) -> Option<&BStr> {
        match self {
            Self::Sam(record) => sam::alignment::Record::name(record),
            Self::Bam(record) => sam::alignment::Record::name(record),
            Self::Cram(record) => sam::alignment::Record::name(record),
        }
    }

    fn flags(&self) -> io::Result<sam::alignment::record::Flags> {
        match self {
            Self::Sam(record) => sam::alignment::Record::flags(record),
            Self::Bam(record) => sam::alignment::Record::flags(record),
            Self::Cram(record) => sam::alignment::Record::flags(record),
        }
    }

    fn reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        match self {
            Self::Sam(record) => sam::alignment::Record::reference_sequence_id(record, header),
            Self::Bam(record) => sam::alignment::Record::reference_sequence_id(record, header),
            Self::Cram(record) => sam::alignment::Record::reference_sequence_id(record, header),
        }
    }

    fn alignment_start(&self) -> Option<io::Result<Position>> {
        match self {
            Self::Sam(record) => sam::alignment::Record::alignment_start(record),
            Self::Bam(record) => sam::alignment::Record::alignment_start(record),
            Self::Cram(record) => sam::alignment::Record::alignment_start(record),
        }
    }

    fn mapping_quality(&self) -> Option<io::Result<sam::alignment::record::MappingQuality>> {
        match self {
            Self::Sam(record) => sam::alignment::Record::mapping_quality(record),
            Self::Bam(record) => sam::alignment::Record::mapping_quality(record),
            Self::Cram(record) => sam::alignment::Record::mapping_quality(record),
        }
    }

    fn cigar(&self) -> Box<dyn sam::alignment::record::Cigar + '_> {
        match self {
            Self::Sam(record) => sam::alignment::Record::cigar(record),
            Self::Bam(record) => sam::alignment::Record::cigar(record),
            Self::Cram(record) => sam::alignment::Record::cigar(record),
        }
    }

    fn mate_reference_sequence_id<'r, 'h: 'r>(
        &'r self,
        header: &'h sam::Header,
    ) -> Option<io::Result<usize>> {
        match self {
            Self::Sam(record) => sam::alignment::Record::mate_reference_sequence_id(record, header),
            Self::Bam(record) => sam::alignment::Record::mate_reference_sequence_id(record, header),
            Self::Cram(record) => {
                sam::alignment::Record::mate_reference_sequence_id(record, header)
            }
        }
    }

    fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        match self {
            Self::Sam(record) => sam::alignment::Record::mate_alignment_start(record),
            Self::Bam(record) => sam::alignment::Record::mate_alignment_start(record),
            Self::Cram(record) => sam::alignment::Record::mate_alignment_start(record),
        }
    }

    fn template_length(&self) -> io::Result<i32> {
        match self {
            Self::Sam(record) => sam::alignment::Record::template_length(record),
            Self::Bam(record) => sam::alignment::Record::template_length(record),
            Self::Cram(record) => sam::alignment::Record::template_length(record),
        }
    }

    fn sequence(&self) -> Box<dyn sam::alignment::record::Sequence + '_> {
        match self {
            Self::Sam(record) => sam::alignment::Record::sequence(record),
            Self::Bam(record) => sam::alignment::Record::sequence(record),
            Self::Cram(record) => sam::alignment::Record::sequence(record),
        }
    }

    fn quality_scores(&self) -> Box<dyn sam::alignment::record::QualityScores + '_> {
        match self {
            Self::Sam(record) => sam::alignment::Record::quality_scores(record),
            Self::Bam(record) => sam::alignment::Record::quality_scores(record),
            Self::Cram(record) => sam::alignment::Record::quality_scores(record),
        }
    }

    fn data(&self) -> Box<dyn sam::alignment::record::Data + '_> {
        match self {
            Self::Sam(record) => sam::alignment::Record::data(record),
            Self::Bam(record) => sam::alignment::Record::data(record),
            Self::Cram(record) => sam::alignment::Record::data(record),
        }
    }
}

impl Default for Record {
    fn default() -> Self {
        Self::Sam(sam::Record::default())
    }
}
