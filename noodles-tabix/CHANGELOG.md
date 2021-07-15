# Changelog

## Unreleased

### Fixed

  * writer: Avoid casts that may truncate.

    Fields that convert to `i32` from other integer types now check whether
    they are in range.

## 0.1.0 - 2021-07-14

  * Initial release.
