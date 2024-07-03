# Changelog

## Unreleased

### Changed

  * bed: Move `Record` to `feature::RecordBuf`.

### Removed

  * bed/feature/record_buf: Remove `str::FromStr` and `fmt::Display`.

    Use a BED record reader and writer, respectively, instead.

## 0.15.0 - 2024-06-17

### Changed

  * bed/record/name: Remove validation on parse ([#271]).

    This allows for more relaxed input in memory and is checked on
    serialization instead.

[#271]: https://github.com/zaeleus/noodles/issues/271

### Fixed

  * bed/record: Fix parse error kind returned from an invalid score ([#270]).

[#270]: https://github.com/zaeleus/noodles/pull/270

## 0.14.0 - 2024-06-06

### Changed

  * bed: Move reader (`Reader`) and writer (`Writer`) to `io` module.

### Deprecated

  * bed: Deprecate `bed::Reader` and `bed::Writer`.

    Use `bed::io::Reader` and `bed::io::Writer`, respectively, instead.

## 0.13.0 - 2024-05-08

### Changed

  * bed: Sync dependencies.

## 0.12.0 - 2024-01-25

### Changed

  * bed: Sync dependencies.

## 0.11.0 - 2023-12-14

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.10.0 - 2023-06-15

### Changed

  * bed: Sync dependencies.

## 0.9.0 - 2023-05-18

### Added

  * bed/record/score: Add a const constructor (`Score::new`), getter
    (`Score::get`), and min (`Score::MIN`) and max (`Score::MAX`) constants.

### Fixed

  * bed/record/name: Fix valid name character set.

    This previously disallowed the tilde (`~`) character.

## 0.8.0 - 2023-03-03

### Changed

  * bed: Sync dependencies.

## 0.7.0 - 2023-02-03

### Added

  * bed/record: Implement `std::error::Error::source` for errors.

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.64.0.

## 0.6.0 - 2022-11-18

### Added

  * bed/record/color: Add getters for color components (`Color::red`,
    `Color::green`, and `Color::blue`).

## 0.5.0 - 2022-10-20

### Changed

  * bed: Sync dependencies.

## 0.4.0 - 2022-08-16

### Changed

  * bed: Raise minimum supported Rust version (MSRV) to 1.57.0.

## 0.3.0 - 2022-06-08

### Added

  * bed/record: Add support for BED7+ (`Record<7>`), BED8+ (`Record<8>`), BED9+
    (`Record<9>`), and BED12+ (`Record<12>`).

## 0.2.0 - 2022-03-29

### Changed

  * bed/record: Wrap start and end positions.

    Positions are now 1-based, inclusive.

## 0.1.0 - 2022-02-17

  * bed: Initial release.
