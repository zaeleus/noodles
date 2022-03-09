# Changelog

## Unreleased

### Added

  * sam: Add an alignment reader trait (`AlignmentReader`).

    This is a generalization over SAM-like alignment formats.

  * sam/alignment_record: Add alignment record fields:

      * flags (`AlignmentRecord::flags`),
      * mapping quality (`AlignmentRecord::mapping_quality`),
      * mate position (`AlignmentRecord::mate_alignment_start`),
      * read name (`AlignmentRecord::read_name`),
      * template length (`AlignmentRecord::template_length`), and
      * quality scores (`AlignmentRecord::quality_scores`).

  * sam/reader: Implement `AlignmentReader`.

  * sam/record/cigar: Add `clear` method.

  * sam/record/cigar: Implement `AsMut<Vec<Op>>`.

  * sam/record/quality_scores/score: Implement `Default` + `Ord` +
    `PartialOrd`.

  * sam/record/quality_scores/score: Add `clear` and `push` methods.

  * sam/record/read_name: Implement `AsRef<[u8]>`, `Hash`, `Ord`, `PartialOrd`,
    and `TryFrom<Vec<u8>>`.

  * sam/record/read_name: Add conversions from `Into<Vec<u8>>`
    (`ReadName::try_new`) and to `Vec<u8>`.

  * sam/record/sequence: Implement `DerefMut`.

  * sam/record/sequence: Add conversion from `Vec<u8>`.

  * sam/record/sequence: Add getting (`Sequence::get`) and setting
    (`Sequence::get_mut`) a base by `Position`.

  * sam/record/sequence/base: Implement `TryFrom<u8>`.

  * sam/record/sequence/base: Implement conversion to `u8`.

### Changed

  * sam: Rename `RecordExt` to `AlignmentRecord`.

  * sam/record: Move `flags`, `mapping_quality`, `read_name`,
    `template_length`, and `quality_scores` to the implementation of
    `AlignmentRecord`.

  * sam/record: Errors that store lengths now use `usize`:

      * `sam::record::parser::SequenceLengthMismatch`,
      * `sam::record::parser::QualityScoresLengthMismatch`,
      * `sam::record::builder::BuildError::SequenceLengthMismatch`, and
      * `sam::record::builder::BuildError::QualityScoresLengthMismatch`.

  * sam/record/cigar/op: Change length to a `usize`.

  * sam/record/cigar/op/kind: Rename `SeqMatch` to `SequenceMatch` and
    `SeqMismatch` to `SequenceMismatch`.

  * sam/record/read_name: Disallow "`*`" as a read name.

    This is already treated as missing (`None`).

### Removed

  * sam/record/read_name: Remove `Deref<Target = String>`.

    Use `AsRef<str>` instead.

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
