// RFC 1952 § 2.3.1
pub(crate) const MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

pub(crate) const MTIME_NONE: u32 = 0;

#[non_exhaustive]
pub(crate) enum CompressionMethod {
    Deflate = 8,
}

#[non_exhaustive]
pub(crate) enum OperatingSystem {
    Unknown = 255,
}
