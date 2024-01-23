//! SAM record field.

use std::io;

use lexical_core::FromLexical;
use noodles_core::Position;

use super::{Bounds, Cigar, Data, Name, QualityScores, ReferenceSequenceName, Sequence};
use crate::Header;

const MISSING: &[u8] = b"*";

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: Vec<u8>,
    pub(crate) bounds: Bounds,
}

impl Fields {
    pub fn name(&self) -> Option<Name<'_>> {
        match &self.buf[self.bounds.name_range()] {
            MISSING => None,
            buf => Some(Name::new(buf)),
        }
    }

    pub fn flags(&self) -> io::Result<u16> {
        let src = &self.buf[self.bounds.flags_range()];
        lexical_core::parse(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }

    pub fn reference_sequence_id(&self, header: &Header) -> Option<io::Result<usize>> {
        self.reference_sequence_name()
            .map(|reference_sequence_name| {
                get_reference_sequence_id(header, reference_sequence_name.as_ref())
            })
    }

    pub fn reference_sequence_name(&self) -> Option<ReferenceSequenceName<'_>> {
        match &self.buf[self.bounds.reference_sequence_name_range()] {
            MISSING => None,
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    pub fn alignment_start(&self) -> Option<io::Result<Position>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.alignment_start_range()] {
            MISSING => None,
            buf => Some(parse_position(buf)),
        }
    }

    pub fn mapping_quality(&self) -> Option<io::Result<u8>> {
        const MISSING: &[u8] = b"255";

        match &self.buf[self.bounds.mapping_quality_range()] {
            MISSING => None,
            buf => Some(parse_int(buf)),
        }
    }

    pub fn cigar(&self) -> Cigar<'_> {
        match &self.buf[self.bounds.cigar_range()] {
            MISSING => Cigar::new(b""),
            buf => Cigar::new(buf),
        }
    }

    pub fn mate_reference_sequence_id(&self, header: &Header) -> Option<io::Result<usize>> {
        self.mate_reference_sequence_name()
            .map(|mate_reference_sequence_name| {
                get_reference_sequence_id(header, mate_reference_sequence_name.as_ref())
            })
    }

    pub fn mate_reference_sequence_name(&self) -> Option<ReferenceSequenceName<'_>> {
        const EQ: &[u8] = b"=";

        match &self.buf[self.bounds.mate_reference_sequence_name_range()] {
            MISSING => None,
            EQ => self.reference_sequence_name(),
            buf => Some(ReferenceSequenceName::new(buf)),
        }
    }

    pub fn mate_alignment_start(&self) -> Option<io::Result<Position>> {
        const MISSING: &[u8] = b"0";

        match &self.buf[self.bounds.mate_alignment_start_range()] {
            MISSING => None,
            buf => Some(parse_position(buf)),
        }
    }

    pub fn template_length(&self) -> io::Result<i32> {
        let buf = &self.buf[self.bounds.template_length_range()];
        parse_int(buf)
    }

    pub fn sequence(&self) -> Sequence<'_> {
        let buf = match &self.buf[self.bounds.sequence_range()] {
            MISSING => b"",
            buf => buf,
        };

        Sequence::new(buf)
    }

    pub fn quality_scores(&self) -> QualityScores<'_> {
        let buf = match &self.buf[self.bounds.quality_scores_range()] {
            MISSING => b"",
            buf => buf,
        };

        QualityScores::new(buf)
    }

    pub fn data(&self) -> Data<'_> {
        let buf = &self.buf[self.bounds.data_range()];
        Data::new(buf)
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            buf: Vec::from(b"*4*0255**00**"),
            bounds: Bounds::default(),
        }
    }
}

fn get_reference_sequence_id(header: &Header, reference_sequence_name: &[u8]) -> io::Result<usize> {
    header
        .reference_sequences()
        .get_index_of(reference_sequence_name)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid reference sequence name",
            )
        })
}

fn parse_position(buf: &[u8]) -> io::Result<Position> {
    parse_int::<usize>(buf).and_then(|n| {
        Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn parse_int<N: FromLexical>(buf: &[u8]) -> io::Result<N> {
    lexical_core::parse(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}
