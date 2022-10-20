# Changelog

## 0.9.0 - 2022-10-20

### Added

  * core/position: Add `const` getter (`Position::get`).

## 0.8.0 - 2022-08-16

### Changed

  * core: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.7.0 - 2022-06-08

### Added

  * core/position: Add `MAX` associated constant.

  * core/region: Add `Interval` container.

    This now represents a closed interval and can now be used separate from
    `Region`.

  * core/region/interval: Add display formatter for `(Unbounded, Included)`.

### Changed

  * core/region: Replace `ParseError::InvalidStartPosition` and
    `ParseError::InvalidEndPosition` with `ParseError::InvalidInterval`.

    The specific errors are now part of `interval::ParseError`.

### Removed

  * core/region: Remove the `Interval` type alias.

    Use the `Interval` struct instead.

## 0.6.0 - 2022-03-29

### Added

  * core: Add a 1-based position wrapper (`Position`).

### Changed

  * core/region: Replace `Region` with `region::Mapped`.

    A `Region` can now only represent a span over a single reference sequence.

  * core/region: Change interval to `Position` bounds.

### Fixed

  * core: Remove dependency on noodles-sam in manifest.

## 0.5.0 - 2022-03-02

### Removed

  * core/region: Remove `from_str_reference_sequences`.

    This is to remove the dependency on noodles-sam. Use `FromStr` instead,
    although it is a less strict parser.

## 0.4.0 - 2022-02-17

### Added

  * core: Set minimum supported Rust version (MSRV) to 1.56.0 in package
    manifest.

## 0.3.4 - 2022-01-27

### Fixed

  * core: Sync dependencies.

## 0.3.3 - 2022-01-13

### Fixed

  * core: Sync dependencies.

## 0.3.2 - 2021-12-09

### Fixed

  * core: Sync dependencies.

## 0.3.1 - 2021-11-18

### Fixed

  * core: Sync dependencies.

## 0.3.0 - 2021-11-11

### Changed

  * core: Update to Rust 2021.

## 0.2.3 - 2021-10-16

### Fixed

  * core: Sync dependencies.

## 0.2.2 - 2021-10-01

### Fixed

  * core: Sync dependencies.

## 0.2.1 - 2021-09-23

### Fixed

  * core: Sync dependencies.

## 0.2.0 - 2021-09-19

### Changed

  * core/region/mapped: Use `std::ops::Bound::cloned` to clone interval bounds.

    This was stabilized in Rust 1.55.0.

## 0.1.2 - 2021-07-30

### Fixed

  * core: Sync dependencies.

## 0.1.1 - 2021-07-21

### Fixed

  * core: Fixed documentation link in package manifest ([#31]).

[#31]: https://github.com/zaeleus/noodles/issues/31

## 0.1.0 - 2021-07-14

  * core: Initial release.
