mod itf8;
mod ltf8;
pub(crate) mod vlq;

use std::{
    io::{self, Read},
    mem, num,
};

pub use self::{
    itf8::{read_itf8, read_itf8_as},
    ltf8::read_ltf8_as,
    vlq::{read_sint7, read_sint7_64, read_uint7, read_uint7_64, read_uint7_as},
};

use crate::file_definition::Version;

/// Reads a signed variable-length integer.
///
/// Uses ITF8 for CRAM 2.x/3.x, sint7 (zigzag-encoded) for CRAM 4.0.
/// Used for codec parameters that are genuinely signed (offsets, alphabet values).
pub fn read_signed_int<R>(reader: &mut R, version: Version) -> io::Result<i32>
where
    R: Read,
{
    if version.uses_vlq() {
        read_sint7(reader)
    } else {
        read_itf8(reader)
    }
}

/// Reads a header-level integer (unsigned VLQ reinterpreted as i32).
///
/// In CRAM 4.0, header-level integers (container and slice headers) use unsigned
/// VLQ encoding. Values that are logically negative (like reference_sequence_id
/// = -1 for unmapped, -2 for multi-ref) are encoded as their unsigned 32-bit
/// representation and reinterpreted as i32 via bitcast.
///
/// For CRAM 2.x/3.x, uses ITF8 which handles both positive and negative natively.
pub fn read_header_int<R>(reader: &mut R, version: Version) -> io::Result<i32>
where
    R: Read,
{
    if version.uses_vlq() {
        let n = read_uint7(reader)?;
        Ok(n as i32)
    } else {
        read_itf8(reader)
    }
}

/// Reads an unsigned variable-length integer as i32.
///
/// Uses ITF8 for CRAM 2.x/3.x, uint7 for CRAM 4.0.
pub fn read_unsigned_int<R>(reader: &mut R, version: Version) -> io::Result<i32>
where
    R: Read,
{
    if version.uses_vlq() {
        let n = read_uint7(reader)?;
        i32::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        read_itf8(reader)
    }
}

/// Reads an unsigned variable-length integer, converting to type `N`.
///
/// Uses ITF8 for CRAM 2.x/3.x, uint7 for CRAM 4.0.
/// For the VLQ path, converts directly from `u32` to `N`,
/// avoiding an unnecessary `i32` intermediate that would reject values > `i32::MAX`.
pub fn read_unsigned_int_as<R, N>(reader: &mut R, version: Version) -> io::Result<N>
where
    R: Read,
    N: TryFrom<u32> + TryFrom<i32>,
    <N as TryFrom<u32>>::Error: Into<Box<dyn std::error::Error + Send + Sync>>,
    <N as TryFrom<i32>>::Error: Into<Box<dyn std::error::Error + Send + Sync>>,
{
    if version.uses_vlq() {
        let n = read_uint7(reader)?;
        N::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        let n = read_itf8(reader)?;
        N::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}

/// Reads a variable-length long integer, converting to type `N`.
///
/// Uses LTF8 for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
pub fn read_long_as<R, N>(reader: &mut R, version: Version) -> io::Result<N>
where
    R: Read,
    N: TryFrom<i64, Error = num::TryFromIntError>,
{
    if version.uses_vlq() {
        let n = read_uint7_64(reader)?;
        let n = i64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        n.try_into()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    } else {
        read_ltf8_as(reader)
    }
}

/// Reads a position value (alignment_start, alignment_span).
///
/// Uses ITF8 (widened to i64) for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
/// CRAM 4.0 header positions use unsigned VLQ (not zigzag-encoded).
pub fn read_position<R>(reader: &mut R, version: Version) -> io::Result<i64>
where
    R: Read,
{
    if version.has_64bit_positions() {
        read_uint7_64(reader).and_then(|n| {
            i64::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    } else {
        read_itf8(reader).map(i64::from)
    }
}

pub(crate) fn read_u8<R>(reader: &mut R) -> io::Result<u8>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u8>()];
    reader.read_exact(&mut buf)?;
    Ok(buf[0])
}

fn read_u16_be<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_be_bytes(buf))
}

pub(crate) fn read_u16_le<R>(reader: &mut R) -> io::Result<u16>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u16>()];
    reader.read_exact(&mut buf)?;
    Ok(u16::from_le_bytes(buf))
}

fn read_u24_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf[1..])?;
    Ok(u32::from_be_bytes(buf))
}

pub(crate) fn read_i32_le<R>(reader: &mut R) -> io::Result<i32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<i32>()];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_u32_be<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_be_bytes(buf))
}

pub(crate) fn read_u32_le<R>(reader: &mut R) -> io::Result<u32>
where
    R: Read,
{
    let mut buf = [0; mem::size_of::<u32>()];
    reader.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}
