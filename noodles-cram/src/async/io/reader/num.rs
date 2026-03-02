mod itf8;
mod ltf8;
mod vlq;

use std::num;

use tokio::io::{self, AsyncRead};

pub use self::{
    itf8::read_itf8,
    ltf8::read_ltf8_as,
    vlq::{read_sint7, read_uint7, read_uint7_64, read_uint7_as},
};

use crate::file_definition::Version;

/// Reads a signed variable-length integer.
///
/// Uses ITF8 for CRAM 2.x/3.x, sint7 (zigzag) for CRAM 4.0.
pub async fn read_signed_int<R>(reader: &mut R, version: Version) -> io::Result<i32>
where
    R: AsyncRead + Unpin,
{
    if version.uses_vlq() {
        read_sint7(reader).await
    } else {
        read_itf8(reader).await
    }
}

/// Reads a header-level integer (unsigned VLQ reinterpreted as i32).
///
/// In CRAM 4.0, all header-level integers (container/block/slice headers) use
/// unsigned VLQ encoding. Values that are logically negative (like
/// reference_sequence_id = -1 for unmapped, -2 for multi-ref) are encoded as
/// their unsigned 32-bit representation and must be reinterpreted as i32.
///
/// For CRAM 2.x/3.x, uses ITF8 which handles both positive and negative natively.
pub async fn read_header_int<R>(reader: &mut R, version: Version) -> io::Result<i32>
where
    R: AsyncRead + Unpin,
{
    if version.uses_vlq() {
        let n = read_uint7(reader).await?;
        Ok(n as i32)
    } else {
        read_itf8(reader).await
    }
}

/// Reads an unsigned variable-length integer, converting to type `N`.
///
/// Uses ITF8 for CRAM 2.x/3.x, uint7 for CRAM 4.0.
/// For the VLQ path, converts directly from `u32` to `N`,
/// avoiding an unnecessary `i32` intermediate that would reject values > `i32::MAX`.
pub async fn read_unsigned_int_as<R, N>(reader: &mut R, version: Version) -> io::Result<N>
where
    R: AsyncRead + Unpin,
    N: TryFrom<u32> + TryFrom<i32>,
    <N as TryFrom<u32>>::Error: Into<Box<dyn std::error::Error + Send + Sync>>,
    <N as TryFrom<i32>>::Error: Into<Box<dyn std::error::Error + Send + Sync>>,
{
    if version.uses_vlq() {
        let n = read_uint7(reader).await?;
        N::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        let n = read_itf8(reader).await?;
        N::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

/// Reads a variable-length long integer, converting to type `N`.
///
/// Uses LTF8 for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
pub async fn read_long_as<R, N>(reader: &mut R, version: Version) -> io::Result<N>
where
    R: AsyncRead + Unpin,
    N: TryFrom<i64, Error = num::TryFromIntError>,
{
    if version.uses_vlq() {
        let n = read_uint7_64(reader).await?;
        let n = i64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        read_ltf8_as(reader).await
    }
}

/// Reads a position value (alignment_start, alignment_span).
///
/// Uses ITF8 (widened to i64) for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
/// CRAM 4.0 header positions use unsigned VLQ (not zigzag-encoded).
pub async fn read_position<R>(reader: &mut R, version: Version) -> io::Result<i64>
where
    R: AsyncRead + Unpin,
{
    if version.has_64bit_positions() {
        read_uint7_64(reader).await.and_then(|n| {
            i64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    } else {
        read_itf8(reader).await.map(i64::from)
    }
}
