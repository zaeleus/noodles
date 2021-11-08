# Changelog

## Unreleased

### Added

  * bam/record: Add mutable getters for flags (`Record::flags_mut`) and mapping
    quality (`Record::mapping_quality_mut`).

  * bam/record/data: Add ordered map methods: `get_index`, `insert`, and `len`.

  * bam/record/data/field: Add fallible conversion to `Vec<u8>`.

  * bam/record/data/field/value/subtype: Add conversion to `u8`.

  * bam/record/sequence: Add append base to sequence (`Sequence::push`).

  * bam/record/quality_scores: Add single score reader (`QualityScores::get`).

### Changed

  * Update to Rust 2021.

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

  * writer/record: Avoid casting data field value array length.

    Converting from `usize` to `i32` now checks whether the value is in range.

## 0.6.0 - 2021-10-02

### Changed

  * Bump minor version due to removal in 0.5.2.

## 0.5.2 - 2021-10-01 

This was erroneously released as a patch release. Use 0.6.0 instead.

### Removed

  * record/data: Remove `Reader`.

    This moves the fields iterator up a module to `bam::record::data`.

## 0.5.1 - 2021-09-23

### Fixed

  * Sync dependencies.

## 0.5.0 - 2021-09-19

### Added

  * bai/async/writer: Add shutdown delegate.

### Deprecated

  * record/data/reader: Deprecate `Reader::read_value_type`.

    Use `noodles_bam::reader::record::data::field::read_value` instead.

### Fixed

  * bai/async: Fix writer not finalizing.

  * bam/reader/record/data/field/value: Avoid casting array length.

    Converting from `i32` to `usize` now checks whether the value is in range.

## 0.4.0 - 2021-08-19

### Changed

  * Update to tokio 1.10.0.

  * async: I/O builders are now owned/consuming builders ([#36]).

    This fixes the terminal method not being able to move out of a mutable
    reference.

  * writer/record: Write lengths as `u32`.

    `block_size` and `l_seq` were changed to unsigned integers in
    [samtools/hts-specs@31d6e44887ae3892472c20d06c15e9a763f3c7c0].

[#36]: https://github.com/zaeleus/noodles/pull/36
[samtools/hts-specs@31d6e44887ae3892472c20d06c15e9a763f3c7c0]: https://github.com/samtools/hts-specs/commit/31d6e44887ae3892472c20d06c15e9a763f3c7c0

### Fixed

  * Define features to enable for Docs.rs.

## 0.3.0 - 2021-08-11

### Added

  * async: Add async reader (`bam::AsyncReader`).

  * async: Add async writer (`bam::AsyncWriter`).

  * bai/async: Add async reader (`bai::AsyncReader`).

  * bai/async: Add async writer (`bai::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

## 0.2.1 - 2021-07-30

### Fixed

  * bai/reader: Return I/O errors when failing to read `n_no_coor`.

    This previously ignored all I/O errors but now only catches
    `UnexpectedEof`.

## 0.2.0 - 2021-07-21

### Added

  * bai/index: Implemented `BinningIndex` for `Index`.

  * bai/index: Added `query` method to find chunks that intersect the given
    region.

  * bai/index/reference_sequence: Implemented `BinningIndexReferenceSequence`
    for `ReferenceSequence`.

  * reader: Accept any `BinningIndex` to query.

### Fixed

  * Fixed documentation link in package manifest ([#31]).

  * bai/reader: Avoid casts that may truncate.

    Fields that convert from `u32` to other integer types now check whether
    they are in range.

  * bai/writer: Avoid casts that may truncate.

    Fields that convert to `u32` from other integer types now check whether
    they are in range.

  * writer/record: Fix bin calculation start position.

    This incorrectly used a 1-based position rather than 0-based.

[#31]: https://github.com/zaeleus/noodles/issues/31

### Removed

  * bai/index/reference_sequence: Removed `Metadata`.

    Use `noodles_csi::index::reference_sequence::Metadata` instead.

## 0.1.0 - 2021-07-14

  * Initial release.
