# Changelog

## 0.21.0 - 2022-10-28

### Added

  * sam/header/record/value/map/tag: Add `AsRef<[u8; LENGTH]>` as a supertrait
    of `Standard`.

### Changed

  * sam/header/parser: Move record parser errors to `record::ParseError`.

    This moves `header::parser::ParseError::InvalidHeader`,
    `header::parser::ParseError::InvalidReferenceSequence`,
    `header::parser::ParseError::InvalidReadGroup`, and
    `header::parser::ParseError::InvalidProgram` to
    `header::record::ParseError`.

### Fixed

  * sam/header/record/value/map/header: Fix serialization when group order
    (`GO`) is set ([#122]).

    This would previously duplicate the group order as the subsort order
    (`SO`).

[#122]: https://github.com/zaeleus/noodles/issues/122

## 0.20.0 - 2022-10-20

### Added

  * sam/alignment_reader: Add reader as a type parameter.

  * sam/record/cigar/op/kind: Add whether an operation kind consumes the read
    (`Kind::consumes_read`) and/or reference (`Kind::consumes_reference`).

  * sam/record/mapping_quality: Add `const` getter (`MappingQuality::get`).

  * sam/reader: Implement `From<R> for Reader<R>`.

### Changed

  * sam: Raise minimum supported Rust version (MSRV) to 1.62.0.

## 0.19.0 - 2022-09-29

### Added

  * sam/header/record/value/map/read_group/platform: Add Element Biosciences
    (`ELEMENT`; [samtools/hts-specs@4b884a4]) and Ultima Genomics (`ULTIMA`;
    [samtools/hts-specs@44b4167]).

   * sam/record/quality_scores/score: Add `MIN` and `MAX` associated constants.

   * sam/record/quality_scores/score: Add const getter (`Score::get`).

[samtools/hts-specs@4b884a4]: https://github.com/samtools/hts-specs/commit/4b884a4bbef181a73c7a3e7d03c91b3fd6371c4d
[samtools/hts-specs@44b4167]: https://github.com/samtools/hts-specs/commit/44b4167209aa0d9a89881dacfab8545f204e4f82

### Changed

  * sam/header: SAM header records parsed from maps are now map values
    (`noodles_sam::header::record::value::Map`).

    This rewraps `header::Header`, `Program`, `ReadGroup`, and
    `ReferenceSequence` as map values `Map<I>`, where `I` is a specialized map
    type. A map is required to have an inner `I` type that can have required
    standard fields and can include optional fields, where the key is a
    nonstandard tag (`tag::Other`).

  * sam/header/record/value/map/reference_sequence: Rename
   `ReferenceSequence::len` to `ReferenceSequence::length`.

### Removed

  * sam/header/record: Remove `Value`.

    Records are fully parsed rather than returning a raw value. Use
    `header::Record` instead.

## 0.18.0 - 2022-08-16

### Changed

  * sam: Raise minimum supported Rust version (MSRV) to 1.59.0.

    This fixes `TryFrom<char> for u8` not being available before 1.59.0
    ([#81]).

[#81]: https://github.com/zaeleus/noodles/pull/81

## 0.17.0 - 2022-07-05

### Added

  * sam/reader: Add query iterator (`Reader::query`).

### Changed

  * sam/alignment: Increase the visibility of `sam::alignment::record`.

    This allows access to the `Builder` struct.

  * sam/record/data/field: `ParseError::InvalidValue` no longer wraps a source.

  * sam/record/data/field/value: `Value::from_str_type` returns `io::Result`.

    This now uses the record data field value parser.

## 0.16.0 - 2022-06-08

### Added

  * sam/alignment: Add a base alignment record (`alignment::Record`).

  * sam/async/writer: Add alignment record writer
    (`Writer::write_alignment_record`).

  * sam/lazy: Add lazy record (`lazy::Record`).

    Lazy records are alignment records that are lazily-evaluated. They fields
    are not necessarily valid, but the buffer is guaranteed to be record-like.

  * sam/reader: Add `Reader::read_lazy_record` to read lazy records.

  * sam/record: Add template length wrapper (`TemplateLength`).

### Changed

  * sam: Replace `Record` with `alignment::Record`.

  * sam/alignment_reader: Return `alignment::Record` from `records` iterator.

  * sam/alignment_writer: Require a concrete `alignment::Record`.

  * sam/async/reader: Read record into a `Record` buffer.

  * sam/async/reader: `Reader::read_record` and `Reader::records` require a
    context (`Header`) when reading a record.

  * sam/async/writer: `Writer::write_record` requires a context (`Header`) when
    writing a record.

  * sam/header/reference_sequence: Change length to a `NonZeroUsize`.

  * sam/reader: Read record into a `Record` buffer.

  * sam/reader: `Reader::read_record` and `Reader::records` require a context
    (`Header`) when reading a record.

  * sam/record/cigar: Change conversion from `Vec<Op>` to be
    fallible.

    Replace usages of `From<Vec<Op>>` with `TryFrom<Vec<Op>>`.

  * sam/record/data/field/value: Validate character values.

    This guarantees the character value is `!`..=`~`.

  * sam/record/data/field/value: Add character value wrapper (`Character`).

  * sam/record/data/field/value: Rename `Value::Char` to `Value::Character`.

  * sam/record/data/field/value: Rename `ParseError::InvalidCharValue` to
    `Value::InvalidCharacterValue`.

  * sam/record/data/field/value/ty: Rename `Type::Char` to `Type::Character`.

  * sam/record/quality_scores: `TryFrom<Vec<u8>>` can now be empty.

  * sam/writer: `Writer::write_record` requires a context (`Header`) when
    writing a record.

  * sam/writer/record: Restrict position to [0, 2^31-1].

    This is defined in _Sequence Alignment/Map Format Specification_
    (2021-06-03) § 1.4 "The alignment section: mandatory fields".

### Deprecated

  * sam/record/cigar: Deprecate `Cigar::read_len` and `Cigar::reference_len`.

    Use `Cigar::read_length` and `Cigar::alignment_span`, respectively,
    instead.

  * sam/record/data/field/value: Deprecate `Value::as_char` and
    `Value::is_char`.

    Use `Value::as_character` and `Value::is_character`, respectively, instead.

### Removed

  * sam: Remove `Record`.

    Use `alignment::Record` instead.

  * sam: Remove `AlignmentRecord`.

    These methods are moved to `alignment::Record`.

  * sam/record: Remove `record::Builder`.

    Use `alignment::record::Builder` instead.

## 0.15.0 - 2022-04-14

### Added

  * sam/record/data/field/value: Implement `TryFrom<char>`.

  * sam/record/sequence: Add conversion to `Vec<Base>`.

  * sam/record/quality_scores: Add conversion to `Vec<Score>`.

  * sam/record/quality_scores: Implement `TryFrom<Vec<u8>>`.

### Changed

  * sam/record/data/field/value: Change conversion from `String` to be
    fallible.

    This replaces `From<String>` with `TryFrom<String>`, which validates the
    string value.

## 0.14.0 - 2022-03-29

### Added

  * sam: Add an alignment reader trait (`AlignmentReader`).

    This is a generalization for reading SAM-like alignment formats.

  * sam: Add an alignment writer trait (`AlignmentWriter`).

    This is a generalization for writing SAM-like alignment formats.

  * sam/alignment_record: Add alignment record fields:

      * cigar (`AlignmentRecord::cigar`),
      * data (`AlignmentRecord::data`),
      * flags (`AlignmentRecord::flags`),
      * mapping quality (`AlignmentRecord::mapping_quality`),
      * mate position (`AlignmentRecord::mate_alignment_start`),
      * read name (`AlignmentRecord::read_name`),
      * sequence (`AlignmentRecord::sequence`),
      * template length (`AlignmentRecord::template_length`), and
      * quality scores (`AlignmentRecord::quality_scores`).

  * sam/header/reference_sequence: Add mutable getter for MD5 checksum
    (`ReferenceSequence::md5_checksum_mut`).

  * sam/reader: Implement `AlignmentReader`.

  * sam/record/data/field/value: Add integer types.

    This means that SAM data field integer values are typed, i.e., it adds
    the `Int8` (`c`), `UInt8` (`C`), `Int16` (`s`), `UInt16` (`S`), and
    `UInt32` (`I`) types. Additionally, it renames `Int` to `Int32`.

    Though these aren't actually used in the SAM format, they can transparently
    be used in SAM and shared among the different alignment formats. An integer
    (`i`) is placed in the smallest type when parsed.

  * sam/record/data/field/value: Add conversions from raw types.

  * sam/record/data/field/value/subtype: Implement `TryFrom<u8>`.

  * sam/record/data/field/value/subtype: Implement conversion to `u8`.

  * sam/record/data/field/value/ty: Implement `Hash` and `TryFrom<u8>`.

  * sam/record/data/field/value/ty: Implement conversion to `u8`.

  * sam/record/cigar: Add `clear` method.

  * sam/record/cigar: Implement `AsMut<Vec<Op>>`.

  * sam/record/mapping_quality: Add `new` method.

  * sam/record/quality_scores/score: Implement `Default` + `Ord` +
    `PartialOrd`.

  * sam/record/quality_scores/score: Add `clear` and `push` methods.

  * sam/record/read_name: Implement `AsRef<[u8]>`, `Hash`, `Ord`, `PartialOrd`,
    and `TryFrom<Vec<u8>>`.

  * sam/record/read_name: Add conversions from `Into<Vec<u8>>`
    (`ReadName::try_new`) and to `Vec<u8>`.

  * sam/record/sequence: Implement `AsRef<[Base]>`, `AsMut<Vec<Base>>`
    and `Index`/`IndexMut` with `Position`.

  * sam/record/sequence: Add conversion from `Vec<u8>`.

  * sam/record/sequence: Add getting (`Sequence::get`) and setting
    (`Sequence::get_mut`) a base by `Position`.

  * sam/record/sequence: Add `clear`, `is_empty`, `len`, and `push` methods.

  * sam/record/sequence/base: Implement `TryFrom<u8>`.

  * sam/record/sequence/base: Implement conversion to `u8`.

  * sam/record/sequence/base: Parse case-insensitive.

  * sam/record/quality_scores: Implement `AsRef<[Score]>`,
    `AsMut<Vec<Score>>`, and `Index`/`IndexMut` with `Position`.

  * sam/record/sequence: Add getting (`QualityScores::get`) and setting
    (`QualityScores::get_mut`) a score by `Position`.

  * sam/record/quality_scores: Add `is_empty` and `len` methods.

  * sam/writer: Implement `AlignmentWriter`.

### Changed

  * sam: Rename `RecordExt` to `AlignmentRecord`.

  * sam/alignment_record: Change alignment start and end to a `Position`.

  * sam/alignment_record: Change alignment span to a `usize`.

  * sam/alignment_record: Remove `Result` from alignment end.

  * sam/record: Move `cigar`, `data`, `flags`, `mapping_quality`,
    `read_name`, `sequence`, `template_length`, and `quality_scores` to the
    implementation of `AlignmentRecord`.

  * sam/record: Errors that store lengths now use `usize`:

      * `sam::record::parser::SequenceLengthMismatch`,
      * `sam::record::parser::QualityScoresLengthMismatch`,
      * `sam::record::builder::BuildError::SequenceLengthMismatch`, and
      * `sam::record::builder::BuildError::QualityScoresLengthMismatch`.

  * sam/record: Change position and mate position to `Position`.

  * sam/record/data/field/value: Rename `Int` to `Int32`.

    As the name implies, the largest value is now `i32::MAX`.

  * sam/record/data/field/value/ty: Rename `Int` to `Int32`.

  * sam/record/cigar: Remove missing state.

    `Cigar` no longer parses `*` as missing. This is now an invalid input.

  * sam/record/cigar/op: Change length to a `usize`.

  * sam/record/cigar/op/kind: Rename `SeqMatch` to `SequenceMatch` and
    `SeqMismatch` to `SequenceMismatch`.

  * sam/record/parser: `ParseError::InvalidPosition` and
    `ParseError::InvalidMatePosition` wrap a `std::num::ParseIntError`.

  * sam/record/read_name: Disallow "`*`" as a read name.

    This is already treated as missing (`None`).

  * sam/record/sequence: Remove missing state.

    `Sequence` no longer parses `*` as missing. This is now an invalid input.

  * sam/record/quality_scores: Remove missing state.

    `QualityScores` no longer parses `*` as missing. This is now parsed as a
    single quality score of 9.

### Removed

  * sam/record: Remove `Position`.

    Use `noodles_core::Position` instead.

  * sam/record/builder: Remove `BuildError`.

    This is no longer used.

  * sam/record/data/field/value: Remove `FromStr`.

    Use `Value::from_str_type` instead.

  * sam/record/read_name: Remove `Deref<Target = String>`.

    Use `AsRef<str>` instead.

  * sam/record/sequence: Remove `Deref<Target = Vec<Base>>`.

    Use `AsRef<[Base]>` instead.

  * sam/record/quality_scores: Remove `Deref<Target = [Score]>`.

    Use `AsRef<[Score]>` instead.

## 0.13.0 - 2022-03-02

### Added

  * sam/record/sequence: Increase the visibility of the `base` module.

## 0.12.0 - 2022-02-17

### Added

  * sam: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

  * sam/async/reader: Add common methods to access the underlying reader:
    `get_ref` and `get_mut`.

  * sam/async/writer: Add missing `get_mut` method.

  * sam/record/position: Implement `Hash` for `Position`.

### Deprecated

  * sam/header/builder: Deprecate `Builder::new`.

    Use `sam::Header::builder` instead.

  * sam/header/header/version: Deprecate `ParseError::MissingMajorVersion` and
    `ParseError::MissingMinorVersion`.

    Use `ParseError::Invalid` instead.

## 0.11.0 - 2022-01-27

### Added

  * sam/record/mapping_quality: Implement `Ord` + `PartialOrd`.

  * sam/record/position: Implement `Ord` + `PartialOrd`.

## 0.10.0 - 2022-01-13

### Added

  * sam/header: Add `clear` method to remove all records from the header.

  * sam/record/data/field/tag: Add base modification probabilities (`ML`) and
    base modifications (`MM`) tags.

  * sam/record/mapping_quality: Add constant for raw missing mapping quality
    (`mapping_quality::MISSING`).

  * sam/record/mapping_quality: Add parser.

### Changed

  * sam/record: Mapping quality is now stored as an `Option`.

    Valid mapping qualities are between 0 and 254, inclusive (`Some`). A
    mapping quality of 255 is considered to be missing (`None`).

  * sam/record/parser: Change `ParseError::InvalidMappingQuality` from wrapping
    `num::ParseIntError` to `mapping_quality::ParseError`.

    Use `mapping_quality::ParseError::Parse` for the `num::ParseIntError`.

### Removed

  * sam/record/mapping_quality: Remove `Deref` and `From<u8>` implementations.

    Conversion is now fallible; use `TryFrom<u8>` instead.

## 0.9.0 - 2021-12-09

### Added

  * sam/writer: Add mutable getter (`Writer::get_mut`) for the underlying
    writer.

### Deprecated

  * sam/record/builder: Deprecate `Builder::new`.

    Use `sam::Record::builder` instead.

  * sam/record/flags: Rename "pair" to "segment".

    The following are now deprecated and listed with their replacements:

      * `Flags::PAIRED` => `Flags::SEGMENTED`
      * `Flags::is_paired` => `Flags::is_segmented`
      * `Flags::PROPER_PAIR` => `Flags::PROPERLY_ALIGNED`
      * `Flags::is_proper_pair` => `Flags::is_properly_aligned`
      * `Flags::READ_1` => `Flags::FIRST_SEGMENT`
      * `Flags::is_read_1` => `Flags::is_first_segment`
      * `Flags::READ_2` => `Flags::LAST_SEGMENT`
      * `Flags::is_read_2` => `Flags::is_last_segment`

## 0.8.1 - 2021-12-02

### Fixed

  * sam: Sync dependencies.

## 0.8.0 - 2021-11-18

### Added

  * sam: Add trait for record extensions (`RecordExt`).

    This includes getting the associated reference sequence
    (`RecordExt::reference_sequence`), alignment start
    (`RecordExt::alignment_start`) and end (`RecordExt::alignment_end`)
    positions, alignment span over the reference sequence
    (`RecordExt::alignment_span`), and mate's associated reference sequence
    (`RecordExt::mate_reference_sequence`).

  * sam/reader: Add common methods to access the underlying reader: `get_ref`
    and `get_mut`.

  * sam/record/reference_sequence_name: Implement `Hash` for
    `ReferenceSequenceName`.

### Changed

  * sam/header/reference_sequence: The reference sequence name (`name`) is now
    wrapped as a `header::reference_sequence::Name`. (This is currently the
    same as a `record::ReferenceSequenceName`.)

  * sam/header/reference_sequence: `ParseError::DuplicateReferenceSequenceName`
    now stores a `header::reference_sequence::Name` rather than a `String`.

## 0.7.0 - 2021-11-11

### Changed

  * sam: Update to Rust 2021.

  * sam/header: Change header, program, read group, reference sequence tag
    representations to `[u8; 2]`.

  * sam/record/data: Specialize the inner ordered map (#49).

    This changes `Data` from wrapping `IndexMap` to a specialized ordered map.
    This improves both parse and access times.

    The interface is familiar but may diverge from the normal collections map
    interface, e.g., keys are passed by value because `Tag` is `Copy` and
    insertion only takes a field because the field tag is inherently the key.

  * sam/record/data/field: Merge `ParseError::MissingTag` and
    `ParseError::MissingValue` into `ParseError::Invalid`.

    The colon (`:`) delimiter must exist to parse the tag and value.

  * sam/record/data/field/tag: Prevent construction of `Tag::Other`.

    Use `Tag::try_from::<[u8; 2]>` instead. This prevents invalid tags from
    being created.

  * sam/record/data/field/value: Merge `ParseError::MissingType` and
    `ParseError::MissingValue` into `ParseError::Invalid`.

    The colon (`:`) delimiter must exist to parse the value.

[#49]: https://github.com/zaeleus/noodles/issues/49

## 0.6.0 - 2021-10-16

### Added

  * sam/record: Add mutable getters for position (`Record::position_mut`), mate
    position (`Record::mate_position_mut`), mate reference sequence name
    (`Record::mate_reference_sequence_name_mut`), and reference sequence name
    (`Record::reference_sequence_name_mut`).

  * sam/record/data/field/tag: Add `ParseError` variants for invalid length and
    character.

  * sam/record/data/field/tag: Implement `Copy` for `Tag` ([#48]).

### Changed

  * sam/record/data: Moved `DuplicateTag` to `ParseError` ([#48]).

    Use `ParseError::DuplicateTag` instead of `ParseError::InvalidData(_)`.

  * sam/record/data/field: `Field::tag` returns a copy rather than a reference
    ([#48]).

  * sam/record/data/field/tag: `Tag::Other` stores a `[u8; 2]` rather than
    `String` ([#46]).

[#46]: https://github.com/zaeleus/noodles/issues/46
[#48]: https://github.com/zaeleus/noodles/pull/48

## 0.5.0 - 2021-10-01

### Added

  * sam/record: Add mutable getters for flags (`Record::flags_mut`; [#39]),
    read name (`Record::read_name_mut`), mapping quality
    (`Record::mapping_quality_mut`), and template length
    (`Record::template_length_mut`).

[#39]: https://github.com/zaeleus/noodles/pull/39

## 0.4.0 - 2021-09-23

### Changed

  * sam/async/reader: Handle CRLF newlines and missing final newline.

  * sam/header/record: Require delimiter split when parsing.

    This ensures the delimiter exists when tokenizing. It removes
    `ParseError::MissingValue` for `ParseError::Invalid` and
    `ParseError::MissingTag` for `ParseError::InvalidField`.

  * sam/header/record: Validate field values (`/[ -~]+/`).

  * sam/reader: Handle CRLF newlines and missing final newline.

## 0.3.0 - 2021-09-19

### Added

  * sam/async: Add async reader (`sam::AsyncReader`).

  * sam/async: Add async writer (`sam::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

  * sam/reader: Add method to return the underlying reader
    (`Reader::into_inner`).

  * sam/writer: Add method to return the underlying writer
    (`Writer::into_inner`).

### Changed

  * sam/header/program/builder: Return error from `build`.

    This previously panicked if the ID was not set.

  * sam/header/read_group/builder: Return error from `build`.

    This previously panicked if the ID was not set.

## 0.2.2 - 2021-08-19

### Fixed

  * sam: Sync dependencies.

## 0.2.1 - 2021-08-11

### Fixed

  * sam: Sync dependencies.

## 0.2.0 - 2021-07-30

### Added

  * sam/header/reference_sequence: Add alternative locus (`AH`) parser.

### Changed

  * sam/header/reference_sequence: Parse alternative locus (`AH`) value.

    This is no longer stored as a raw `String`.

## 0.1.1 - 2021-07-21

### Fixed

  * sam: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * sam: Initial release.
