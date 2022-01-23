# Changelog

## Unreleased

### Deprecated

  * fasta/record/definition: Deprecate
    `ParseError::MissingReferenceSequenceName`.

    Use `ParseError::MissingName` instead.

## 0.5.2 - 2022-01-13

### Fixed

  * fasta: Sync dependencies.

## 0.5.1 - 2021-12-09

### Fixed

  * fasta: Sync dependencies.

## 0.5.0 - 2021-12-02

### Added

  * fasta/record/definition: Accept `Into<String>` for name.

### Changed

  * fasta/record: Wrap sequence.

    Use `Sequence::as_ref` to get the underlying list.

## 0.4.1 - 2021-11-18

### Fixed

  * fasta: Sync dependencies.

## 0.4.0 - 2021-11-11

### Changed

  * fasta: Update to Rust 2021.

## 0.3.1 - 2021-10-16

### Fixed

  * fasta: Sync dependencies.

## 0.3.0 - 2021-10-01

### Deprecated

  * fasta/fai/record: `Record::reference_sequence_name` is now `Record::name`.

    FASTA records are not necessarily reference sequences.

## 0.2.4 - 2021-09-23

### Fixed

  * fasta: Sync dependencies.

## 0.2.3 - 2021-09-19

### Fixed

  * fasta: Sync dependencies.

## 0.2.2 - 2021-08-19

### Fixed

  * fasta: Sync dependencies.

## 0.2.1 - 2021-08-11

### Fixed

  * fasta: Sync dependencies.

## 0.2.0 - 2021-07-30

### Changed

  * fasta/record: Rename `reference_sequence_name` to `name`.

    FASTA records are not necessarily reference sequences.

## 0.1.1 - 2021-07-21

### Fixed

  * fasta: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * fasta: Initial release.
