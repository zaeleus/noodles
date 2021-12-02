# Changelog

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
