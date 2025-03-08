# Changelog

## 0.17.0 - 2025-03-08

### Changed

  * core: Raise minimum supported Rust version (MSRV) to 1.81.0.

## 0.16.0 - 2025-01-19

### Changed

  * core: Raise minimum supported Rust version (MSRV) to 1.73.0.

## 0.15.0 - 2024-05-08

### Changed

  * core/region: Wrap name as a `BString`.

    This is nearly identical to its previous `Vec<u8>` type but allows for
    better ergonomics.

### Removed

  * core: Remove `Error`.

    This is unused.

## 0.14.0 - 2024-01-25

### Changed

  * core/region: Change name to a byte string (`Vec<u8>`).

## 0.13.0 - 2023-12-14

### Changed

  * core: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.12.0 - 2023-06-15

### Added

  * core/region/interval: Add position inclusion query (`Interval::contains`).

## 0.11.0 - 2023-03-03

### Added

  * core/position: Make `Position::checked_add` `const`.

## 0.10.0 - 2023-02-03

### Added

  * core: Add a generalized error type (`Error`).

  * core/region: Implement `std::error::Error::source` for errors.

### Changed

  * core: Raise minimum supported Rust version (MSRV) to 1.64.0.

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
