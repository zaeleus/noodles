# Changelog

## 0.84.0 - 2025-11-13

### Changed

  * bam/io/reader/query: Change `Query` to a reader ([#265]).

    `Query` is now a reader adapter, rather than an iterator. Use
    `Query::records` for an iterator.

[#265]: https://github.com/zaeleus/noodles/pull/265

## 0.83.0 - 2025-08-25

### Changed

  * bam: Sync dependencies.

## 0.82.0 - 2025-07-12

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.81.0 - 2025-05-29

### Added

  * bam/record/quality_scores: Add `QualityScores::iter` ([#338]).

    This implementation is also provided by
    `sam::alignment::record::QualityScores`.

[#338]: https://github.com/zaeleus/noodles/pull/338

### Changed

  * bam/record/fields: Return empty quality scores when missing.

    This previously returned a representation of missing quality scores.

### Removed

  * bam: Remove deprecated items.

    The following items are removed:

      * `AsyncReader` (deprecated since 0.68.0; use `r#async::io::Reader`
        instead),
      * `AsyncWriter` (0.68.0; `r#async::io::Writer`),
      * `bai::AsyncReader` (0.68.0; `bai::r#async::io::Reader`),
      * `bai::AsyncWriter` (0.68.0; `bai::r#async::io::Writer`),
      * `bai::Reader` (0.68.0; `bai::io::Reader`),
      * `bai::Writer` (0.68.0; `bai::io::Writer`),
      * `bai::r#async::Reader` (0.68.0; `bai::r#async::io::Reader`),
      * `bai::r#async::Writer` (0.68.0; `bai::r#async::io::Writer`),
      * `bai::r#async::read` (0.73.0; `bai::r#async::fs::read`),
      * `bai::r#async::write` (0.73.0; `bai::r#async::fs::write`),
      * `bai::read` (0.73.0; `bai::fs::read`), and
      * `bai::write` (0.73.0; `bai::fs::write`).

## 0.80.0 - 2025-05-16

### Changed

  * bam: Sync dependencies.

## 0.79.0 - 2025-04-13

### Changed

  * bam: Sync dependencies.

## 0.78.0 - 2025-04-06

### Changed

  * bam/record/codec/encoder/data/field/value: Validate string and hex values.

    String values must now match `[ !-~]*`; and hex values,
    `([0-9A-F][0-9A-F])*`.

## 0.77.0 - 2025-03-08

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.81.0.

### Fixed

  * bam/record/sequence/iter: Fix size hint calculation.

## 0.76.0 - 2025-02-06

### Changed

  * bam: Sync dependencies.

## 0.75.0 - 2025-01-24

### Added

  * bam/fs: Add indexer ([#320]).

[#320]: https://github.com/zaeleus/noodles/issues/320

## 0.74.0 - 2025-01-23

### Changed

  * bam: Sync dependencies.

## 0.73.0 - 2025-01-19

### Added

  * bam/io/reader: Add a header reader adapter (`header::Reader`).

    This can be used to read the raw fields of a BAM header, e.g., the raw SAM
    header.

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.73.0.

  * bam/bai: Move convenience functions (`read` and `write`) to `fs` module.

### Deprecated

  * bam/bai: Deprecate `read` and `write`.

    Use `bai::fs::read` and `bai::fs::write`, respectively, instead.

### Fixed

  * bam/record/codec/encoder/quality_scores: Replace usage of unstable std
    function ([#315]).

    `std::iter::repeat_n` was stabilized in Rust 1.82.0, but the current MSRV
    is set to 1.70.0.

[#315]: https://github.com/zaeleus/noodles/issues/315

## 0.72.0 - 2024-12-20

### Changed

  * bam: Sync dependencies.

## 0.71.0 - 2024-12-12

### Changed

  * bam/io/reader: Use records iterator instead of record buffers iterator for
    alignment records ([#313]).

  * bam/record/sequence: Increase the visibility of `Iter`.

    `Sequence::iter` now returns `Iter` instead of an implementation of
    `Iterator.`

[#313]: https://github.com/zaeleus/noodles/discussions/313

## 0.70.0 - 2024-11-07

### Changed

  * bam: Sync dependencies.

### Fixed

  * bam/record/codec/encoder/sequence: Encode base as case-insensitive.

## 0.69.0 - 2024-10-22

### Changed

  * bam/record/codec/decoder/reference_sequence_id: Remove validation against
    the reference sequence dictionary.

    Validation of the reference sequence ID still occurs on look up.

## 0.68.0 - 2024-09-26

### Added

  * bam/async/io/reader: Add stream for querying unmapped records
    (`Reader::query_unmapped`).

### Changed

  * bam/bai: Move reader (`Reader`) and writer (`Writer`) to `io` module.

  * bam/bai/async: Move reader (`Reader`) and writer (`Writer`) to `io` module.

### Deprecated

  * bam: Deprecate async re-exports (`AsyncReader` and `AsyncWriter`).

    Use `bam::r#async::io::Reader` and `bam::r#async::io::Writer`
    instead.

  * bam/bai: Deprecate `Reader` and `Writer`.

    Use `bai::io::Reader` and `bai::io::Writer`, respectively, instead.

  * bam/bai/async: Deprecate `Reader` and `Writer`.

    Use `bai::r#async::io::Reader` and `bai::r#async::io::Writer`,
    respectively, instead.

## 0.67.0 - 2024-09-04

### Added

  * bam/io/writer/builder: Add build from writer
    (`Builder::build_from_writer`).

### Changed

  * bam/record/codec/decoder: Remove panic when reading bin on EOF.

    This now returns a `DecodeError` if the input reaches EOF.

## 0.66.0 - 2024-08-04

### Added

  * bam/record/sequence: Add get base at index (`Sequence::get`) and checked
    split at (`Sequence::split_at_checked`).

## 0.65.0 - 2024-07-14

### Changed

  * bam: Update to bit-vec 0.7.0.

## 0.64.0 - 2024-06-17

### Added

  * bam/bai: Add common methods to access the underlying I/O.

## 0.63.0 - 2024-05-16

### Changed

  * bam: Sync dependencies.

## 0.62.0 - 2024-05-08

### Changed

  * bam: Sync dependencies.

## 0.61.0 - 2024-05-02

### Removed

  * bam/io/reader: Remove `Reader::virtual_position` and `Reader::seek`.

    Use the inner reader instead.

## 0.60.0 - 2024-04-22

### Changed

  * bam: Sync dependencies.

### Fixed

  * bam/record/name: Trim trailing `NUL` when converting to an alignment record
    name buffer ([#254]).

[#254]: https://github.com/zaeleus/noodles/issues/254

## 0.59.0 - 2024-04-04

### Changed

  * bam: Sync dependencies.

## 0.58.0 - 2024-03-28

### Changed

  * bam/io/reader/record: Disallow partial block sizes.

    A partial block size previously would return EOF. This is now more strict,
    requiring a complete block size or EOF.

## 0.57.0 - 2024-03-12

### Changed

  * bam: Sync dependencies.

### Fixed

  * bam/record/data/field/value/array/values: Fix value count calculation
    ([#244]).

[#244]: https://github.com/zaeleus/noodles/issues/244

## 0.56.0 - 2024-02-15

### Changed

  * bam: Sync dependencies.

## 0.55.0 - 2024-02-08

### Changed

  * bam: Sync dependencies.

## 0.54.1 - 2024-02-04

### Fixed

  * bam/record/fields: Fix consuming CIGAR operations buffer when there are 2
    operations ([#230]).

[#230]: https://github.com/zaeleus/noodles/issues/230

## 0.54.0 - 2024-02-01

### Changed

  * bam/record/codec/encoder: Allow unsupported bin IDs to overflow ([#229]).

[#229]: https://github.com/zaeleus/noodles/issues/229

## 0.53.0 - 2024-01-25

### Added

  * bam/record: Add flags (`record::Flags`), mapping quality
    (`record::MappingQuality`), position (`record::Position`), reference
    sequence ID (`record::ReferenceSequenceId`), and template length
    (`record::TemplateLength`) wrappers.

  * bam/record: Implement `sam::alignment::Record`.

  * bam/record/data/field/value/array: Add values wrapper (`Values`).

### Changed

  * bam: Move lazy record to record.

    Use `bam::Record` instead of `bam::lazy::Record`.

  * bam: Move readers (`Reader` and `IndexedReader`) and writer (`Writer`) to
    `io` module.

  * bam/async/io: Consider binary reference sequences as part of the
    header.

    When reading, if the SAM header has a reference sequence dictionary, it
    must match the binary reference sequences; otherwise, the binary reference
    sequences are added to the SAM header.

    When writing, binary reference sequences are written from the SAM header.

  * bam/async/io/reader: Change `Reader::read_header` to return a parsed header
    (`sam::Header`).

    This no longer returns a raw string.

  * bam/io/reader: Rename "record" to "record buf" and "lazy record" to
    "record".

    This changes the following:

      * `Reader::read_record` => `Reader::read_record_buf`,
      * `Reader::records` => `Reader::record_bufs`,
      * `Records` => `RecordBufs`,
      * `Reader::read_lazy_record` => `Reader::read_record`,
      * `Reader::lazy_records` => `Reader::records`, and
      * `LazyRecords` => `Records`.

  * bam/io/reader: Return an iterator over `bam::Record` for queries
    (`Reader::query` and `Reader::query_unmapped`).

    `Reader::query_unmapped` no longer needs a reference to the header.

  * bam/record: Rename `ReadName` to `Name`.

    This also changes `Record::read_name` to `Record::name`.

  * bam/record/data/field: Replace `Value` with
    `sam::alignment::record::data::field::Value`.

  * bam/record/data/field/value: Replace `Array` with
    `sam::alignment::record::data::field::value::Array`.

### Fixed

  * bam/record: Discard skip length when matching overflowing CIGAR.

  * bam/record/codec/encoder: Use alignment span for `m` in overflowing CIGAR.

### Removed

  * bam/async/io: Remove `Reader::read_reference_sequences` and
    `Writer::write_reference_sequences`.

    These are now considered as part of the header when calling reading or
    writing.

## 0.52.0 - 2023-12-14

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.70.0.

  * bam/async/reader: Accept `csi::BinningIndex` for querying.

  * bam/bai: Define `bai::Index` as `binning_index::Index<LinearIndex>`.

  * bam/reader: Accept `csi::BinningIndex` for querying.

## 0.51.0 - 2023-11-14

### Added

  * bam/bai: Add `Index` type alias for `csi::Index`.

### Removed

  * bam/bai/reader: Remove `read_header`.

    This is now read when reading the index (`Reader::read_index`).

  * bam/bai/writer: Remove `write_header`.

    This is now written when writing the index (`Writer::write_index`).

## 0.50.0 - 2023-11-13

### Changed

  * bam: Sync dependencies.

## 0.49.1 - 2023-11-02

### Fixed

  * bam/reader/header: Handle optional `NUL` padding in SAM header text.

## 0.49.0 - 2023-10-26

### Changed

  * bam: Sync dependencies.

## 0.48.0 - 2023-10-19

### Changed

  * bam: Sync dependencies.

## 0.47.0 - 2023-10-12

### Added

  * bam/lazy/record/read_name: Add `ReadName::as_bytes` to return the buffer
    without the trailing `NUL` terminator.

  * bam/lazy/record/sequence: Add bases iterator (`Sequence::iter`).

## 0.46.0 - 2023-09-21

### Changed

  * bam: Sync dependencies.

## 0.45.0 - 2023-09-14

### Changed

  * bam/bai/async/reader: Disallow duplicate bin IDs.

  * bam/lazy/record: Change `lazy::Record::flags` to be infallible.

## 0.44.0 - 2023-08-31

### Added

  * bam/lazy/record: Increase the visibility of the `data` module.

  * bam/lazy/record/data: Add lookup by tag (`Data::get`).

  * bam/lazy/record/data/field/value: Add `Value::ty` to return the value
    type.

### Changed

  * bam/lazy/record: Change `lazy::Record::mapping_quality` to be infallible.

  * bam/reader: Require a SAM header when querying for unmapped records
    (`Reader::query_unmapped`).

    It's possible for a chunk to include mapped records, which are subsequently
    filtered out, but they do require the associated header to decode.

### Removed

  * bam/reader: Remove `UnmappedRecords` iterator.

    This is replaced by an `impl Trait`. Use `impl Iterator<Item =
    io::Result<Record>>` instead.

## 0.43.0 - 2023-08-24

### Added

  * bam/lazy/record: Add read name wrapper (`lazy::record::ReadName`).

  * bam/lazy/record/cigar: Add op iterator (`lazy::record::Cigar::iter`).

  * bam/lazy/record/data: Add field iterator (`lazy::record::Data::iter`).

### Changed

  * bam/lazy: Increase the visibility of `lazy::record`.

  * bam/record/codec/decoder/data/field/value/array: Read length as `u32`.

    This used to be defined as an `int32_t` but was changed to a `uint32_t` in
    [samtools/hts-specs@31d6e44].

[samtools/hts-specs@31d6e44]: https://github.com/samtools/hts-specs/commit/31d6e44887ae3892472c20d06c15e9a763f3c7c0

## 0.42.0 - 2023-08-17

### Changed

  * bam: Sync dependencies.

## 0.41.0 - 2023-08-03

### Changed

  * bam/reader/header: Parse header line by line.

    The header parser can now build a `sam::Header` by parsing a raw header
    line by line. This makes it so that it is no longer required to read the
    entire raw header into memory before parsing.

## 0.40.0 - 2023-07-27

### Changed

  * bam: Sync dependencies.

## 0.39.0 - 2023-07-20

### Changed

  * bam: Sync dependencies.

## 0.38.0 - 2023-07-06

### Changed

  * bam: Sync dependencies.

## 0.37.0 - 2023-06-29

### Changed

  * bam/bai/reader: Disallow duplicate bin IDs.

## 0.36.0 - 2023-06-15

### Changed

  * bam: Sync dependencies.

### Fixed

  * bam/record/codec/encoder/sequence: Fix writing an empty sequence when the
    read length is nonempty ([#176]).

[#176]: https://github.com/zaeleus/noodles/issues/176

## 0.35.0 - 2023-06-08

### Removed

  * bam/reader: Reduce the visibility of `record`.

    `bam::reader::record::data::field::get_value` is not intended to be part of
    the public API.

  * bam/writer: Reduce the visibility of `record`.

    `bam::writer::record::data::field::put_value` is not intended to be part of
    the public API.

## 0.34.0 - 2023-06-01

### Changed

  * bam: Sync dependencies.

## 0.33.0 - 2023-05-18

### Changed

  * bam: Sync dependencies.

## 0.32.0 - 2023-05-11

### Changed

  * bam: Sync dependencies.

## 0.31.0 - 2023-05-04

### Changed

  * bam/reader/record/data: Disallow duplicate tags.

## 0.30.0 - 2023-04-27

### Added

  * bam/writer: Add builder (`writer::Builder`).

### Fixed

  * bam/reader/record/position: Fix reading position at max position (2^31-1).

    This would previously overflow and error instead of returning a properly
    normalized `Some(Position)`.

## 0.29.0 - 2023-04-06

### Added

  * bam/indexed_reader: Add getter for index (`IndexedReader::index`).

  * bam/indexed_reader/builder: Add `Builder::build_from_reader`.

  * bam/reader: Add builder (`bam::reader::Builder`).

### Changed

  * bam/bai: Replace `Index` with `noodles_csi::Index`.

  * bam/indexed_reader/builder: Attempt to load either an associated BAM index
    (`<src>.bai`) or CSI (`<src>.csi`), in that order.

  * bam/reader: Change `Reader::read_header` to return a parsed header
    (`sam::Header`).

    This no longer returns a raw string.

  * bam/reader: Consider binary reference sequences as part of the header.

    If the SAM header has a reference sequence dictionary, it must match the
    binary reference sequences; otherwise, the binary reference sequences are
    added to the SAM header.

  * bam/writer: Consider binary reference sequences as part of the header.

    The binary reference sequences was now written as part of the header. They
    are taken from the reference sequence dictionary of the given header.

### Removed

  * bam/reader: Remove `Reader::read_reference_sequences`.

    This is now considered as part of the header when calling
    `Reader::read_header`.

  * bam/writer: Remove `Writer::write_reference_sequences`.

    This is now considered as part of the header when calling
    `Writer::write_header`.

## 0.28.0 - 2023-03-14

### Changed

  * bam: Sync dependencies.

## 0.27.0 - 2023-03-03

### Added

  * bam/reader/record: Resolve CIGAR with `CG` data field value when `CIGAR` is
    a flag.

    When `CIGAR` is `kSmN` (`k` = `l_seq`, `m` = associated reference sequence
    length), the decoder will attempt to resolve the actual CIGAR from the `CG`
    data value. See § 4.2.2 "`N_CIGAR_OP` field" (2022-08-22).

  * bam/writer/record: Write CIGAR to `CG` data field when it has more than
    65535 operations.

    This sets the CIGAR field as a flag for a decoder to read the actual CIGAR
    from the `CG` data field instead. See § 4.2.2 "`N_CIGAR_OP` field"
    (2022-08-22).

### Changed

  * bam/async/reader: Change `Reader::query` to receive a header
    (`sam::Header`) rather than a reference sequence dictionary
    (`sam::header::ReferenceSequences`).

  * bam/async/reader: Require `Reader::read_record` and `Reader::records` to
    receive a header (`sam::Header`).

  * bam/reader: Change `Reader::query` to receive a header (`sam::Header`)
    rather than a reference sequence dictionary
    (`sam::header::ReferenceSequences`).

  * bam/reader: Require `Reader::read_record` and `Reader::records` to receive
    a header (`sam::Header`).

  * bam/reader/record: Validate reference sequence IDs.

    This includes both the reference sequence ID (`ref_id`) and mate reference
    sequence ID (`next_ref_id`) fields. Reference sequence IDs must be either
    missing (`-1`) or less than the number of reference sequence dictionary
    entries (`n_ref`). See § 4.2 "The BAM format" (2022-08-22).

## 0.26.0 - 2023-02-03

### Changed

  * bam: Raise minimum supported Rust version (MSRV) to 1.64.0.

### Removed

  * bam/bai/index: Remove `Index::unplaced_unmapped_read_count`.

    This was deprecated in noodles-bam 0.2.0. Use
    `Index::unplaced_unmapped_record_count` instead.

  * bam/record: Remove `reference_sequence_id`.

    This was deprecated in noodles-bam 0.19.0. Use `-1` instead of
    `reference_sequence_id::UNMAPPED`.

## 0.25.1 - 2022-11-29

### Changed

  * bam: Sync dependencies.

## 0.25.0 - 2022-11-18

### Added

  * bam/writer/record/sequence: Add read length validation.

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
