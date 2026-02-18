# Changelog

## 0.9.0 - 2026-02-18

### Changed

  * refget: Raise minimum supported Rust version (MSRV) to 1.88.0.

  * refget: Update to reqwest 0.13.1.

    noodles-refget no longer enables a TLS backend. Select one by adding
    reqwest as a direct dependency.

## 0.8.0 - 2025-07-12

### Changed

  * refget: Raise minimum supported Rust version (MSRV) to 1.85.0.

## 0.7.0 - 2025-03-08

### Changed

  * refget: Raise minimum supported Rust version (MSRV) to 1.81.0.

  * refget/sequence: Handle response failures ([#322]).

    Response errors now wrap the HTTP client error as `Error::Response`.

[#322]: https://github.com/zaeleus/noodles/issues/332

## 0.6.0 - 2025-01-19

### Changed

  * refget: Raise minimum supported Rust version (MSRV) to 1.73.0.

## 0.5.0 - 2024-05-08

### Changed

  * refget: Sync dependencies.

## 0.4.0 - 2024-04-04

### Changed

  * refget: Update to reqwest 0.12.2.

## 0.3.0 - 2024-01-25

### Changed

  * refget: Increase the visibility of the `sequence` module.

## 0.2.0 - 2023-12-14

### Changed

  * refget: Raise minimum supported Rust version (MSRV) to 1.70.0.

## 0.1.0 - 2023-08-03

  * refget: Initial release.
