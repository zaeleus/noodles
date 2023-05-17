use noodles_core::Position;
use noodles_sam as sam;

use crate::record::{Features, Flags, NextMateFlags};

#[derive(Default)]
pub struct Chunk {
    ids: Vec<u64>,
    bam_bit_flags: Vec<sam::record::Flags>,
    cram_bit_flags: Vec<Flags>,
    reference_sequence_ids: Vec<Option<usize>>,
    read_lengths: Vec<usize>,
    alignment_starts: Vec<Option<Position>>,
    read_group_ids: Vec<Option<usize>>,
    read_names: Vec<Option<sam::record::ReadName>>,
    next_mate_bit_flags: Vec<NextMateFlags>,
    next_fragment_reference_sequence_ids: Vec<Option<usize>>,
    next_mate_alignment_starts: Vec<Option<Position>>,
    template_sizes: Vec<i32>,
    distances_to_next_fragment: Vec<Option<usize>>,
    tags: Vec<sam::record::Data>,
    bases: Vec<sam::record::Sequence>,
    features: Vec<Features>,
    mapping_qualities: Vec<Option<sam::record::MappingQuality>>,
    quality_scores: Vec<sam::record::QualityScores>,
}

impl Chunk {
    pub fn resize(&mut self, new_len: usize) {
        self.ids.resize(new_len, 0);
        self.bam_bit_flags
            .resize(new_len, sam::record::Flags::default());
        self.cram_bit_flags.resize(new_len, Flags::default());
        self.reference_sequence_ids.resize(new_len, None);
        self.read_lengths.resize(new_len, 0);
        self.alignment_starts.resize(new_len, None);
        self.read_group_ids.resize(new_len, None);
        self.read_names.resize(new_len, None);
        self.next_mate_bit_flags
            .resize(new_len, NextMateFlags::default());
        self.next_fragment_reference_sequence_ids
            .resize(new_len, None);
        self.next_mate_alignment_starts.resize(new_len, None);
        self.template_sizes.resize(new_len, 0);
        self.distances_to_next_fragment.resize(new_len, None);
        self.tags.resize(new_len, sam::record::Data::default());
        self.bases.resize(new_len, sam::record::Sequence::default());
        self.features.resize(new_len, Features::default());
        self.mapping_qualities.resize(new_len, None);
        self.quality_scores
            .resize(new_len, sam::record::QualityScores::default());
    }

    #[rustfmt::skip]
    pub fn records_mut(&mut self) -> impl Iterator<Item = RecordMut<'_>> {
        self.ids
            .iter_mut()
            .zip(&mut self.bam_bit_flags)
            .zip(&mut self.cram_bit_flags)
            .zip(&mut self.reference_sequence_ids)
            .zip(&mut self.read_lengths)
            .zip(&mut self.alignment_starts)
            .zip(&mut self.read_group_ids)
            .zip(&mut self.read_names)
            .zip(&mut self.next_mate_bit_flags)
            .zip(&mut self.next_fragment_reference_sequence_ids)
            .zip(&mut self.next_mate_alignment_starts)
            .zip(&mut self.template_sizes)
            .zip(&mut self.distances_to_next_fragment)
            .zip(&mut self.tags)
            .zip(&mut self.bases)
            .zip(&mut self.features)
            .zip(&mut self.mapping_qualities)
            .zip(&mut self.quality_scores)
            .map(
                |(((((((((((((((((id, bam_bit_flags), cram_bit_flags), reference_sequence_id), read_length), alignment_start), read_group_id), read_name), next_mate_bit_flags), next_fragment_reference_sequence_id), next_mate_alignment_start), template_size), distance_to_next_fragment), tags), bases), features), mapping_quality), quality_scores)| {
                    RecordMut {
                        id,
                        bam_bit_flags,
                        cram_bit_flags,
                        reference_sequence_id,
                        read_length,
                        alignment_start,
                        read_group_id,
                        read_name,
                        next_mate_bit_flags,
                        next_fragment_reference_sequence_id,
                        next_mate_alignment_start,
                        template_size,
                        distance_to_next_fragment,
                        tags,
                        bases,
                        features,
                        mapping_quality,
                        quality_scores,
                    }
                },
            )
    }

    #[rustfmt::skip]
    pub fn into_record_bufs(self) -> impl Iterator<Item = crate::Record> {
        self.ids
            .into_iter()
            .zip(self.bam_bit_flags.into_iter())
            .zip(self.cram_bit_flags.into_iter())
            .zip(self.reference_sequence_ids.into_iter())
            .zip(self.read_lengths.into_iter())
            .zip(self.alignment_starts.into_iter())
            .zip(self.read_group_ids.into_iter())
            .zip(self.read_names.into_iter())
            .zip(self.next_mate_bit_flags.into_iter())
            .zip(self.next_fragment_reference_sequence_ids.into_iter())
            .zip(self.next_mate_alignment_starts.into_iter())
            .zip(self.template_sizes.into_iter())
            .zip(self.distances_to_next_fragment.into_iter())
            .zip(self.tags.into_iter())
            .zip(self.bases.into_iter())
            .zip(self.features.into_iter())
            .zip(self.mapping_qualities.into_iter())
            .zip(self.quality_scores.into_iter())
            .map(
                |(((((((((((((((((id, bam_bit_flags), cram_bit_flags), reference_sequence_id), read_length), alignment_start), read_group_id), read_name), next_mate_bit_flags), next_fragment_reference_sequence_id), next_mate_alignment_start), template_size), distance_to_next_fragment), tags), bases), features), mapping_quality), quality_scores)| {
                    crate::Record {
                        id,
                        bam_bit_flags,
                        cram_bit_flags,
                        reference_sequence_id,
                        read_length,
                        alignment_start,
                        read_group: read_group_id,
                        read_name,
                        next_mate_bit_flags,
                        next_fragment_reference_sequence_id,
                        next_mate_alignment_start,
                        template_size,
                        distance_to_next_fragment,
                        tags,
                        bases,
                        features,
                        mapping_quality,
                        quality_scores,
                    }
                },
            )
    }
}

pub struct RecordMut<'a> {
    pub id: &'a mut u64,
    pub bam_bit_flags: &'a mut sam::record::Flags,
    pub cram_bit_flags: &'a mut Flags,
    pub reference_sequence_id: &'a mut Option<usize>,
    pub read_length: &'a mut usize,
    pub alignment_start: &'a mut Option<Position>,
    pub read_group_id: &'a mut Option<usize>,
    pub read_name: &'a mut Option<sam::record::ReadName>,
    pub next_mate_bit_flags: &'a mut NextMateFlags,
    pub next_fragment_reference_sequence_id: &'a mut Option<usize>,
    pub next_mate_alignment_start: &'a mut Option<Position>,
    pub template_size: &'a mut i32,
    pub distance_to_next_fragment: &'a mut Option<usize>,
    pub tags: &'a mut sam::record::Data,
    pub bases: &'a mut sam::record::Sequence,
    pub features: &'a mut Features,
    pub mapping_quality: &'a mut Option<sam::record::MappingQuality>,
    pub quality_scores: &'a mut sam::record::QualityScores,
}