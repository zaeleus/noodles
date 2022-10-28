# Changelog

## 0.24.0 - 2022-10-28

### Changed

  * bam: Sync dependencies.

## 0.23.0 - 2022-10-20

### Added

  * bam: Add an indexed reader (`IndexedReader`).

  * bam/lazy/record: Add conversion to `sam::alignment::Record`.

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.62.0.

### Fixed

  * bam/bai/index/builder: Ensure reference sequence IDs are increasing when
    adding records.

## 0.22.0 - 2022-09-29

### Added

  * bam/lazy/record: Implement `AsRef<[u8]>` and `TryFrom<Vec<u8>>` (#113).

[#113]: https://github.com/zaeleus/noodles/issues/113

## 0.21.0 - 2022-08-16

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.59.0.

### Removed

  * bam/async/reader: Remove `Builder` and `Reader::builder`.

    Use `bgzf::async::reader::Builder` and construct with
    `bam::async::Reader::from` instead.

  * bam/async/writer: Remove `Builder` and `Writer::builder`.

    Use `bgzf::async::writer::Builder` and construct with
    `bam::async::Writer::from` instead.

## 0.20.0 - 2022-07-05

### Added

  * bam/lazy/record: Add data, CIGAR, quality scores, and
    sequence wrappers.

    These allow the raw data to be accessed. Use conversion methods (`TryFrom`)
    on these wrappers to parse the data.

### Changed

  * bam/bai/index/reference_sequence/bin: Change bin ID to a `usize`.

## 0.19.0 - 2022-06-08

### Added

  * bam/async/reader: Add `Reader::read_lazy_record` and `Reader::lazy_records`
    to read lazy records.

  * bam/lazy: Add lazy record (`lazy::Record`).

    Lazy records are alignment records that are lazily-evaluated. They fields
    are not necessarily valid, but the buffer is guaranteed to be record-like.

  * bam/reader: Add `Reader::read_lazy_record` and `Reader::lazy_records` to
    read lazy records.

  * bam/reader/record: Validate the reference sequence ID.

    This is done by checking whether an entry exists in the reference sequence
    dictionary.

### Changed

  * bam/async/writer: `Writer::write_record` requires a context (`sam::Header`)
    when writing a record.

  * bam: Replace `Record` with `sam::alignment::Record`.

  * bam/reader/query: `Query` no longer has a generic type parameter for the
    query interval.

  * bam/reader/record/data/field/value: Validate character values.

  * bam/writer: `Writer::write_record` requires a context (`sam::Header`) when
    writing a record.

### Deprecated

  * bam/record: Deprecate the `reference_sequence_id` module.

    This includes the `reference_sequence_id::UNMAPPED` constant.

### Removed

  * bam: Remove `Record`.

    Use `sam::alignment::Record` instead.

  * bam/async/writer: Remove `Writer::write_sam_record`.

    Use `Writer::write_alignment_record` instead.

  * bam/record: Remove `record::Builder`.

    Use `sam::alignment::record::Builder` instead.

  * bam/writer: Remove `Writer::write_sam_record`.

    Use `Writer::write_alignment_record` instead.

## 0.18.0 - 2022-04-14

### Added

  * bam/async/writer: Add alignment record writer
    (`Writer::write_alignment_record`).

  * bam/record: Add mutable getter for template length
    (`Record::template_length_mut`).

### Deprecated

  * bam/async/writer: Deprecate `Writer::write_sam_record`.

    Use `Writer::write_alignment_record` instead.

## 0.17.0 - 2022-03-29

### Added

  * bam/reader: Implement `sam::AlignmentReader`.

  * bam/record: Add mutable getter for position (`Record::position_mut`) and
    mate position (`Record::mate_position`).

  * bam/writer: Implement `sam::AlignmentWriter`.

### Changed

  * bam/bai/index/reference_sequence: `ReferenceSequence::query` returns an
    `io::Error` instead of `QueryError`.

  * bam/record: Move `cigar`, `data`, `flags`, `mapping_quality`,
    `read_name`, `sequence`, `template_length`, and `quality_scores` to the
    implementation of `sam::AlignmentRecord`.

  * bam/record: Wrap read name as `sam::record::ReadName`.

    Read names are now guaranteed to be valid. The read name "*" is now
    considered missing (`None`).

  * bam/record: Replace `QualityScores` with `sam::record::QualityScores`.

  * bam/record: Replace `Sequence` with `sam::record::Sequence`.

  * bam/record: Replace `Cigar` with `sam::record::Cigar`.

  * bam/record: Replace `Data` with `sam::record::Data`.

  * bam/record: Change position and mate position to `Position`.

  * bam/record: Change reference sequence ID and mate reference sequence ID to
    an `Option<usize>`.

  * bam/record/data/field: Replace `Value` with
    `sam::record::data::field::Value`.

  * bam/reader/record: Clear quality scores if all scores are missing.

    This previously filled the quality scores with scores of 255, but this just
    signals that the quality scores are missing.

  * bam/reader/record/data/field/value: Rename `read_value` to `get_value`.

    This now takes a buffer instead of a reader.

### Deprecated

  * bam/writer: Deprecate `Writer::write_sam_record`.

    Use `Writer::write_alignment_record` instead.

### Fixed

  * bam/async/writer/record: Calculate bin number for record.

  * bam/bai/index/builder: Handle indexing only unplaced, unmapped records.

  * bam/bai/index/reference_sequence: Ensure the start position is not out of
    range for a query (`2^29 - 1`).

  * bam/bai/index/reference_sequence/builder: Calculate bin number for record.

  * bam/record/builder: Fix default read name to have no NUL terminator.

  * bam/writer/record: Calculate bin number for record.

    If the start position or CIGAR is changed, the bin number in the record may
    be incorrect.

### Removed

  * bam/record: Remove bin number (`Record::bin`).

    This can be calculated from the alignment start and end instead.

  * bam/record: Remove `ReferenceSequenceId`.

    This was deprecated in noodles-bam 0.16.0. Use `usize` instead.

  * bam/record: Remove `QualityScores`.

    Use `sam::record::QualityScores` instead.

  * bam/record: Remove `Sequence`.

    Use `sam::record::Sequence` instead.

  * bam/record: Remove `Cigar`.

    Use `sam::record::Cigar` instead.

  * bam/record: Remove `Data`.

    Use `sam::record::Data` instead.

  * bam/record/builder: Remove `BuildError`.

    This is no longer used.

  * bam/record/data/field: Remove `Value`, `value::Type`, and `value::Subtype`.

    Use `sam::record::data::field::Value`,
    `sam::record::data::field::value::Type`, and
    `sam::record::data::field::value::Subtype`, respectively, instead.

## 0.16.0 - 2022-03-02

### Added

  * bam/async/writer: Add common methods to access the underlying writer:
    `get_ref`, `get_mut`, and `into_inner`.

  * bam/record/sequence: Add fallible conversion from `&[u8]`
    (`TryFrom<&[u8]>`) ([#76]).

[#76]: https://github.com/zaeleus/noodles/issues/76

### Changed

  * bam/record/reference_sequence_id: Change underlying type to a `usize`.

    Use `From<usize>` instead of `From<i32>`.

### Deprecated

  * bam/record/reference_sequence_id: Deprecate conversion from `i32`.

    Convert the raw reference sequence ID to a `usize` first, and use
    `ReferenceSequenceId::from` instead.

## 0.15.0 - 2022-02-17

### Added

  * bam: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

  * bam/record: Add mutable getter for read name (`Record::read_name_mut`).

  * bam/record/reference_sequence_id: Implement `Hash + Ord + PartialOrd` for
    `ReferenceSequenceId`.

  * bam/record/cigar: Add `clear` method to remove all CIGAR operations.

### Changed

  * bam/reader/record: Ensure `l_read_name` is > 0.

    The raw read name will always be at least length 1 because the NUL
    terminator must exist.

  * bam/record: `Record::read_name` now returns a byte string sans the NUL
    terminator.

    This is no longer treated as a `CStr`.

## 0.14.0 - 2022-01-27

### Changed

  * bam/async/reader: Parse record from buffer.

    This previously used an async reader when parsing the record but now just
    reads from the block buffer.

## 0.13.0 - 2022-01-13

### Added

  * bam/async/reader: Add conversion from `R` to `Reader<R>`.

  * bam/async/reader: Add common methods to access the underlying reader:
    `get_ref`, `get_mut`, and `into_inner`.

  * bam/async/writer: Add conversion from `W` to `Writer<W>`.

  * bam/record: Mapping quality is now stored as an `Option`.

    Valid mapping qualities are between 0 and 254, inclusive (`Some`). A
    mapping quality of 255 is considered to be missing (`None`).

## 0.12.0 - 2021-12-16

### Added

  * bam/record/data: Add conversion from `Vec<Field>` to `Data`.

### Fixed

  * bam/async/reader: Remove stray debug output when reading reference
    sequences.

## 0.11.0 - 2021-12-09

### Added

  * bam/reader: Add conversion from `R` to `Reader<R>`.

  * bam/reader: Add common methods to access the underlying reader: `get_ref`,
    `get_mut`, and `into_inner`.

  * bam/writer: Add conversion from `W` to `Writer<W>`.

  * bam/writer: Add common methods to access the underlying writer: `get_mut`
    and `into_inner`.

## 0.10.0 - 2021-12-02

### Added

  * bam/record: Add builder (`bam::record::Builder`).

  * bam/record: Add mutable getters for reference sequence ID
    (`Record::reference_sequence_id_mut`) and mate reference sequence ID
    (`Record::mate_reference_sequence_id_mut`).

  * bam/record/cigar: Add conversion from `Vec<Op>` to `Cigar`.

  * bam/record/cigar/op: Add conversion from `Op` to `sam::record::cigar::Op`.

  * bam/record/data: Add method to retrieve a field by tag (`Data::get`).

  * bam/record/sequence: Add conversion from `Vec<Base>` to `Sequence`.

  * bam/record/quality_scores: Add conversion from `Vec<Score>` to
    `QualityScores`.

### Fixed

  * bam: Require tokio's `fs` feature as a dependency ([#62]).

  * bam/async/writer/record: Fix `n_cigar_op` value ([#58]).

  * bam/async/writer/record: Fix encoding SEQ and QUAL fields ([#59]).

  * bam/writer/record: Fix `n_cigar_op` value ([#58]).

    This used an old calculation when `Cigar` stored a raw byte buffer.

  * bam/writer/record: Fix encoding SEQ and QUAL fields ([#59]).

    This validates the lengths to ensure they're equal or only `QUAL` is
    missing.

  * bam/writer/sam_record: Return error when the sequence is empty and quality
    scores exist.

[#58]: https://github.com/zaeleus/noodles/issues/58
[#59]: https://github.com/zaeleus/noodles/issues/59
[#62]: https://github.com/zaeleus/noodles/issues/62

## 0.9.0 - 2021-11-18

### Added

  * bam/record: Implement `sam::RecordExt`.

## 0.8.0 - 2021-11-11

### Added

  * bam/record: Add mutable getters for flags (`Record::flags_mut`) and mapping
    quality (`Record::mapping_quality_mut`).

  * bam/record/data: Add ordered map methods: `get_index`, `insert`, `keys`,
    and `len`.

  * bam/record/data/field: Add fallible conversion to `Vec<u8>`.

  * bam/record/data/field/value/subtype: Add conversion to `u8`.

  * bam/record/sequence: Add append base to sequence (`Sequence::push`).

  * bam/record/quality_scores: Add single score reader (`QualityScores::get`).

### Changed

  * bam: Update to Rust 2021.

  * bam/record: `bam::Record` is no longer backed by a contiguous buffer.

    Fields are read individually when reading the record. `bam::Record` no
    longer implements `Deref<Target = [u8]>`. `Cigar`, `Sequence`, and
    `QualityScores`, `Data` now own their data.

  * bam/async/writer/builder: The compression level wrapper changed from
    `flate2::Compression` to `noodles_bgzf::writer::CompressionLevel`.

### Deprecated

  * bam/record/cigar: Deprecate `Cigar::new`.

    Use `Cigar::from::<Vec<u32>>` instead.

  * bam/record/cigar/op: Deprecate `TryFrom<&[u8]>`.

    Use `TryFrom<u32>` instead.

  * bam/record/data: Deprecate `Data::fields`.

    Use `Data::values` instead.

  * bam/record/sequence: Deprecate `Sequence::base_count`.

    Use `Sequence::len` instead.

  * bam/record/quality_scores: Deprecate `QualityScores::new`.

    Use `QualityScores::from::<Vec<u8>>` instead.

### Removed

  * bam/record/data: Remove `Data::new`.

    Use `Data::try_from::<Vec<u8>>` instead.

## 0.7.0 - 2021-10-16

### Changed

  * bam/record/data: Moved `DuplicateTag` to `ParseError` ([#48]).

    Use `ParseError::DuplicateTag` instead of `ParseError::InvalidData(_)`.

  * bam/record/data/field: `Field::tag` returns a copy rather than a reference
    ([#48]).

[#48]: https://github.com/zaeleus/noodles/pull/48

### Fixed

  * bam/writer/record: Avoid casting data field value array length.

    Converting from `usize` to `i32` now checks whether the value is in range.

## 0.6.0 - 2021-10-02

### Changed

  * bam: Bump minor version due to removal in 0.5.2.

## 0.5.2 - 2021-10-01 

This was erroneously released as a patch release. Use 0.6.0 instead.

### Removed

  * bam/record/data: Remove `Reader`.

    This moves the fields iterator up a module to `bam::record::data`.

## 0.5.1 - 2021-09-23

### Fixed

  * bam: Sync dependencies.

## 0.5.0 - 2021-09-19

### Added

  * bam/bai/async/writer: Add shutdown delegate.

### Deprecated

  * bam/record/data/reader: Deprecate `Reader::read_value_type`.

    Use `noodles_bam::reader::record::data::field::read_value` instead.

### Fixed

  * bam/bai/async: Fix writer not finalizing.

  * bam/reader/record/data/field/value: Avoid casting array length.

    Converting from `i32` to `usize` now checks whether the value is in range.

## 0.4.0 - 2021-08-19

### Changed

  * bam: Update to tokio 1.10.0.

  * bam/async: I/O builders are now owned/consuming builders ([#36]).

    This fixes the terminal method not being able to move out of a mutable
    reference.

  * bam/writer/record: Write lengths as `u32`.

    `block_size` and `l_seq` were changed to unsigned integers in
    [samtools/hts-specs@31d6e44887ae3892472c20d06c15e9a763f3c7c0].

[#36]: https://github.com/zaeleus/noodles/pull/36
[samtools/hts-specs@31d6e44887ae3892472c20d06c15e9a763f3c7c0]: https://github.com/samtools/hts-specs/commit/31d6e44887ae3892472c20d06c15e9a763f3c7c0

### Fixed

  * bam: Define features to enable for Docs.rs.

## 0.3.0 - 2021-08-11

### Added

  * bam/async: Add async reader (`bam::AsyncReader`).

  * bam/async: Add async writer (`bam::AsyncWriter`).

  * bam/bai/async: Add async reader (`bai::AsyncReader`).

  * bam/bai/async: Add async writer (`bai::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.1 - 2021-07-30

### Fixed

  * bam/bai/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * bam/bai/index: Implemented `BinningIndex` for `Index`.

  * bam/bai/index: Added `query` method to find chunks that intersect the given
    region.

  * bam/bai/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

  * bam/reader: Accept any `BinningIndex` to query.

### Fixed

  * bam: Fixed documentation link in package manifest ([#31]).

  * bam/bai/reader: Avoid casts that may truncate.

    Fields that convert from `u32` to other integer types now check whether
    they are in range.

  * bam/bai/writer: Avoid casts that may truncate.

    Fields that convert to `u32` from other integer types now check whether
    they are in range.

  * bam/writer/record: Fix bin calculation start position.

    This incorrectly used a 1-based position rather than 0-based.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * bam/bai/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * bam: Initial release.
