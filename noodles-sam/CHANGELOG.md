# Changelog

## Unreleased

### Changed

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

  * async/reader: Handle CRLF newlines and missing final newline.

  * header/record: Require delimiter split when parsing.

    This ensures the delimiter exists when tokenizing. It removes
    `ParseError::MissingValue` for `ParseError::Invalid` and
    `ParseError::MissingTag` for `ParseError::InvalidField`.

  * header/record: Validate field values (`/[ -~]+/`).

  * reader: Handle CRLF newlines and missing final newline.

## 0.3.0 - 2021-09-19

### Added

  * async: Add async reader (`sam::AsyncReader`).

  * async: Add async writer (`sam::AsyncWriter`).

    Async I/O can be enabled with the `async` feature.

  * reader: Add method to return the underlying reader (`Reader::into_inner`).

  * writer: Add method to return the underlying writer (`Writer::into_inner`).

### Changed

  * header/program/builder: Return error from `build`.

    This previously panicked if the ID was not set.

  * header/read_group/builder: Return error from `build`.

    This previously panicked if the ID was not set.

## 0.2.2 - 2021-08-19

### Fixed

  * Sync dependencies.

## 0.2.1 - 2021-08-11

### Fixed

  * Sync dependencies.

## 0.2.0 - 2021-07-30

### Added

  * header/reference_sequence: Add alternative locus (`AH`) parser.

### Changed

  * header/reference_sequence: Parse alternative locus (`AH`) value.

    This is no longer stored as a raw `String`.

## 0.1.1 - 2021-07-21

### Fixed

  * Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * Initial release.
