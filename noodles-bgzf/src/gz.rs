// RFC 1952 § 2.3.1
pub(crate) const MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

pub(crate) const MTIME_NONE: u32 = 0;

// ID1 (1) + ID2 (1) + CM (1) + FLG (1) + MTIME (4) + XLF (1) + OS (1)
pub(crate) const HEADER_SIZE: usize = 10;

// CRC32 (4) + ISIZE (4)
pub(crate) const TRAILER_SIZE: usize = 8;

#[non_exhaustive]
pub(crate) enum CompressionMethod {
    Deflate = 8,
}

#[non_exhaustive]
pub(crate) enum OperatingSystem {
    Unknown = 255,
}
