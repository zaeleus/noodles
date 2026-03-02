use std::fmt::Write;

use noodles_core::Position;

use super::Feature;
use crate::container::compression_header::preservation_map::{
    SubstitutionMatrix, substitution_matrix::Base,
};

/// Computed MD string and NM (edit distance) values.
pub struct MdNm {
    /// The MD tag string (mismatched and deleted reference bases).
    pub md: String,
    /// The NM tag value (edit distance to the reference).
    pub nm: u32,
}

/// Computes the MD string and NM edit distance from CRAM features and a reference sequence.
///
/// `features` should be sorted by position.
/// `reference_sequence` is the full reference contig (1-based indexing via alignment_start).
/// `alignment_start` is the 1-based start position of the alignment on the reference.
/// `read_length` is the length of the read.
/// `substitution_matrix` is used to resolve Substitution feature codes to actual bases.
pub fn compute_md_nm(
    features: &[Feature<'_>],
    reference_sequence: &[u8],
    alignment_start: Position,
    read_length: usize,
    substitution_matrix: &SubstitutionMatrix,
) -> MdNm {
    let ref_start = usize::from(alignment_start);
    let mut md = String::new();
    let mut nm: u32 = 0;
    let mut match_count: usize = 0;

    // We walk through the alignment using read position (1-based) and ref offset.
    // Features are sorted by read position. Between features, bases match.
    let mut read_pos: usize = 1; // 1-based position in the read
    let mut ref_offset: usize = 0; // offset from alignment_start on the reference

    for feature in features {
        let feat_pos = usize::from(feature.position()); // 1-based in read

        match feature {
            // Companion features (Scores, QualityScore) can share a position with a
            // preceding multi-base feature that already advanced read_pos past that position.
            // HardClip and Padding also don't affect position tracking.
            Feature::Scores { .. }
            | Feature::QualityScore { .. }
            | Feature::HardClip { .. }
            | Feature::Padding { .. } => continue,
            _ => {
                debug_assert!(
                    feat_pos >= read_pos,
                    "feature position {feat_pos} is before current read position {read_pos}; features must be sorted by position"
                );
            }
        }

        match feature {
            Feature::Substitution { code, .. } => {
                // Bases between last position and this one are matches
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                // Get the reference base
                let ref_base_u8 = reference_sequence
                    .get(ref_start - 1 + ref_offset)
                    .copied()
                    .unwrap_or(b'N');
                let ref_base = Base::try_from(ref_base_u8).unwrap_or(Base::N);
                let _read_base = substitution_matrix.get(ref_base, *code);

                // Flush match count, emit ref base
                write!(md, "{match_count}").unwrap();
                match_count = 0;
                md.push(ref_base_u8.to_ascii_uppercase() as char);
                nm += 1;

                read_pos = feat_pos + 1;
                ref_offset += 1;
            }
            Feature::ReadBase { base, .. } => {
                // ReadBase stores a single base + quality; it represents a match/mismatch
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                let ref_base_u8 = reference_sequence
                    .get(ref_start - 1 + ref_offset)
                    .copied()
                    .unwrap_or(b'N');

                if !base.eq_ignore_ascii_case(&ref_base_u8) {
                    // Mismatch
                    write!(md, "{match_count}").unwrap();
                    match_count = 0;
                    md.push(ref_base_u8.to_ascii_uppercase() as char);
                    nm += 1;
                } else {
                    match_count += 1;
                }

                read_pos = feat_pos + 1;
                ref_offset += 1;
            }
            Feature::Bases { bases, .. } => {
                // Stretch of bases (match or mismatch)
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                for (i, &read_base) in bases.iter().enumerate() {
                    let ref_base_u8 = reference_sequence
                        .get(ref_start - 1 + ref_offset + i)
                        .copied()
                        .unwrap_or(b'N');

                    if !read_base.eq_ignore_ascii_case(&ref_base_u8) {
                        write!(md, "{match_count}").unwrap();
                        match_count = 0;
                        md.push(ref_base_u8.to_ascii_uppercase() as char);
                        nm += 1;
                    } else {
                        match_count += 1;
                    }
                }

                read_pos = feat_pos + bases.len();
                ref_offset += bases.len();
            }
            Feature::Insertion { bases, .. } => {
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                nm += bases.len() as u32;
                read_pos = feat_pos + bases.len();
                // Insertion does not consume reference
            }
            Feature::InsertBase { .. } => {
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                nm += 1;
                read_pos = feat_pos + 1;
                // InsertBase does not consume reference
            }
            Feature::Deletion { len, .. } => {
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                write!(md, "{match_count}").unwrap();
                match_count = 0;
                md.push('^');

                for i in 0..*len {
                    let ref_base_u8 = reference_sequence
                        .get(ref_start - 1 + ref_offset + i)
                        .copied()
                        .unwrap_or(b'N');
                    md.push(ref_base_u8.to_ascii_uppercase() as char);
                }

                nm += *len as u32;
                ref_offset += len;
                read_pos = feat_pos;
                // Deletion does not consume read
            }
            Feature::ReferenceSkip { len, .. } => {
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                // N in CIGAR — skip reference bases, no MD entry
                ref_offset += len;
                read_pos = feat_pos;
            }
            Feature::SoftClip { bases, .. } => {
                // matches_before accounts for aligned bases preceding the clip,
                // advancing ref_offset. The clip itself consumes read but not reference.
                let matches_before = feat_pos - read_pos;
                match_count += matches_before;
                ref_offset += matches_before;

                read_pos = feat_pos + bases.len();
            }
            // Scores, QualityScore, HardClip, Padding are handled by `continue` above.
            Feature::Scores { .. }
            | Feature::QualityScore { .. }
            | Feature::HardClip { .. }
            | Feature::Padding { .. } => unreachable!(),
        }
    }

    // Remaining bases after last feature are matches
    if read_pos <= read_length {
        match_count += read_length - read_pos + 1;
    }

    // Flush final match count
    write!(md, "{match_count}").unwrap();

    MdNm { md, nm }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use super::*;

    #[test]
    fn test_compute_md_nm_all_match() {
        let features = [];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            4,
            &SubstitutionMatrix::default(),
        );
        assert_eq!(result.md, "4");
        assert_eq!(result.nm, 0);
    }

    #[test]
    fn test_compute_md_nm_with_deletion() {
        let features = [Feature::Deletion {
            position: Position::try_from(3).unwrap(),
            len: 2,
        }];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            4,
            &SubstitutionMatrix::default(),
        );
        assert_eq!(result.md, "2^GT2");
        assert_eq!(result.nm, 2);
    }

    #[test]
    fn test_compute_md_nm_with_insertion() {
        let features = [Feature::Insertion {
            position: Position::try_from(3).unwrap(),
            bases: Cow::Borrowed(b"TT"),
        }];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        // Read: A C TT G T (6 bases), alignment span = 4
        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            6,
            &SubstitutionMatrix::default(),
        );
        assert_eq!(result.md, "4");
        assert_eq!(result.nm, 2);
    }

    #[test]
    fn test_compute_md_nm_with_read_base_mismatch() {
        let features = [Feature::ReadBase {
            position: Position::try_from(2).unwrap(),
            base: b'T',
            quality_score: 30,
        }];
        let reference = b"ACGT";
        let alignment_start = Position::try_from(1).unwrap();

        // Position 2 is C in ref, T in read -> mismatch
        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            4,
            &SubstitutionMatrix::default(),
        );
        assert_eq!(result.md, "1C2");
        assert_eq!(result.nm, 1);
    }

    #[test]
    fn test_compute_md_nm_with_substitution() {
        // Default matrix: for ref base C (index 1), code 0b10 => T
        let sub_matrix = SubstitutionMatrix::default();
        let code = sub_matrix.find(Base::C, Base::T);

        let features = [Feature::Substitution {
            position: Position::try_from(2).unwrap(),
            code,
        }];
        let reference = b"ACGT";
        let alignment_start = Position::try_from(1).unwrap();

        // Position 2 is C in ref, substituted to T -> mismatch
        let result = compute_md_nm(&features, reference, alignment_start, 4, &sub_matrix);
        assert_eq!(result.md, "1C2");
        assert_eq!(result.nm, 1);
    }

    #[test]
    fn test_compute_md_nm_combined_features() {
        // CIGAR: 2M 1I 1M 2D 2M = read length 6, alignment span 7
        // Read positions: 1=A 2=C 3=X(ins) 4=G 5=A 6=C
        // Ref positions:  1=A 2=C         3=G 4=T 5=A 6=C 7=G
        // Deletion of ref bases T,A (positions 4-5) happens at read position 5
        // (after the 1M match at read position 4, before the final 2M).
        let sub_matrix = SubstitutionMatrix::default();

        let features = [
            Feature::Insertion {
                position: Position::try_from(3).unwrap(),
                bases: Cow::Borrowed(b"X"),
            },
            Feature::Deletion {
                position: Position::try_from(5).unwrap(),
                len: 2,
            },
        ];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(&features, reference, alignment_start, 6, &sub_matrix);
        // 2 matches, insertion (NM+1), 1 match, deletion of TA (NM+2), 2 matches
        assert_eq!(result.md, "3^TA2");
        assert_eq!(result.nm, 3);
    }

    #[test]
    fn test_compute_md_nm_with_soft_clip() {
        // CIGAR: 2S 3M = read length 5, alignment span 3
        // Soft clip at read start: positions 1-2 are clipped
        let features = [Feature::SoftClip {
            position: Position::try_from(1).unwrap(),
            bases: Cow::Borrowed(b"NN"),
        }];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            5,
            &SubstitutionMatrix::default(),
        );
        // Soft clips don't affect MD/NM, 3 remaining bases match
        assert_eq!(result.md, "3");
        assert_eq!(result.nm, 0);
    }

    #[test]
    fn test_compute_md_nm_companion_features() {
        // Verify that companion features (Scores after Bases) don't trigger debug_assert
        let features = [
            Feature::Bases {
                position: Position::try_from(1).unwrap(),
                bases: Cow::Borrowed(b"AC"),
            },
            Feature::Scores {
                position: Position::try_from(1).unwrap(),
                quality_scores: Cow::Borrowed(&[40, 35]),
            },
        ];
        let reference = b"ACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            4,
            &SubstitutionMatrix::default(),
        );
        assert_eq!(result.md, "4");
        assert_eq!(result.nm, 0);
    }

    #[test]
    fn test_compute_md_nm_with_reference_skip() {
        // CIGAR: 2M 3N 2M = read length 4, alignment span 7
        // Reference skip (N in CIGAR) skips reference bases without MD/NM impact.
        let features = [Feature::ReferenceSkip {
            position: Position::try_from(3).unwrap(),
            len: 3,
        }];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            4,
            &SubstitutionMatrix::default(),
        );
        // 2 matches, skip 3 ref bases, 2 matches
        assert_eq!(result.md, "4");
        assert_eq!(result.nm, 0);
    }

    #[test]
    fn test_compute_md_nm_with_insert_base() {
        // CIGAR: 2M 1I 2M = read length 5, alignment span 4
        // InsertBase is a single-base insertion (vs multi-base Insertion).
        let features = [Feature::InsertBase {
            position: Position::try_from(3).unwrap(),
            base: b'X',
        }];
        let reference = b"ACGTACGT";
        let alignment_start = Position::try_from(1).unwrap();

        let result = compute_md_nm(
            &features,
            reference,
            alignment_start,
            5,
            &SubstitutionMatrix::default(),
        );
        // 2 matches, insertion (NM+1), 2 matches
        assert_eq!(result.md, "4");
        assert_eq!(result.nm, 1);
    }
}
