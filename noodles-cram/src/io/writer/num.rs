mod itf8;
mod ltf8;
mod vlq;

use std::io::{self, Write};

pub use self::{
    itf8::write_itf8,
    ltf8::write_ltf8,
    vlq::{uint7_size_of, write_sint7, write_sint7_64, write_uint7, write_uint7_64},
};

use crate::file_definition::Version;

/// Writes an unsigned variable-length integer (as i32).
///
/// Uses ITF8 for CRAM 2.x/3.x, uint7 for CRAM 4.0.
pub fn write_int<W>(writer: &mut W, version: Version, value: i32) -> io::Result<()>
where
    W: Write,
{
    if version.uses_vlq() {
        let n = u32::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(writer, n)
    } else {
        write_itf8(writer, value)
    }
}

/// Writes a header-level integer (unsigned VLQ from i32).
///
/// In CRAM 4.0, header-level integers use unsigned VLQ encoding. Values that
/// are logically negative (like -1 for unmapped, -2 for multi-ref) are
/// reinterpreted as their unsigned 32-bit representation via bitcast.
///
/// For CRAM 2.x/3.x, uses ITF8 which handles both positive and negative natively.
pub fn write_header_int<W>(writer: &mut W, version: Version, value: i32) -> io::Result<()>
where
    W: Write,
{
    if version.uses_vlq() {
        write_uint7(writer, value as u32)
    } else {
        write_itf8(writer, value)
    }
}

/// Writes a signed variable-length integer.
///
/// Uses ITF8 for CRAM 2.x/3.x, sint7 (zigzag-encoded) for CRAM 4.0.
pub fn write_signed_int<W>(writer: &mut W, version: Version, value: i32) -> io::Result<()>
where
    W: Write,
{
    if version.uses_vlq() {
        write_sint7(writer, value)
    } else {
        write_itf8(writer, value)
    }
}

/// Writes a variable-length long integer.
///
/// Uses LTF8 for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
pub fn write_long<W>(writer: &mut W, version: Version, value: i64) -> io::Result<()>
where
    W: Write,
{
    if version.uses_vlq() {
        let n = u64::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7_64(writer, n)
    } else {
        write_ltf8(writer, value)
    }
}

/// Writes a position value (alignment_start, alignment_span).
///
/// Uses ITF8 (narrowed from i64) for CRAM 2.x/3.x, uint7_64 for CRAM 4.0.
pub fn write_position<W>(writer: &mut W, version: Version, value: i64) -> io::Result<()>
where
    W: Write,
{
    if version.has_64bit_positions() {
        let n = u64::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7_64(writer, n)
    } else {
        let n = i32::try_from(value).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_itf8(writer, n)
    }
}

/// Returns the encoded size of a variable-length integer for the given version.
///
/// Paired with `write_int` for sizing content block fields (`content_id`,
/// `compressed_size`, `uncompressed_size`), which are always non-negative.
/// For CRAM 4.0 (VLQ), `n` is bitcast to `u32` to compute the size.
/// Negative values produce incorrect sizes and `write_int` will reject them.
pub fn int_size_of(version: Version, n: i32) -> usize {
    if version.uses_vlq() {
        uint7_size_of(n as u32)
    } else {
        itf8_size_of(n)
    }
}

/// Returns the encoded size of an ITF8 integer in bytes.
pub fn itf8_size_of(n: i32) -> usize {
    if n >> (8 - 1) == 0 {
        1
    } else if n >> (16 - 2) == 0 {
        2
    } else if n >> (24 - 3) == 0 {
        3
    } else if n >> (32 - 4) == 0 {
        4
    } else {
        5
    }
}

pub(crate) fn write_u8<W>(writer: &mut W, n: u8) -> io::Result<()>
where
    W: Write,
{
    writer.write_all(&[n])
}

pub(crate) fn write_u16_le<W>(writer: &mut W, n: u16) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

pub(crate) fn write_i32_le<W>(writer: &mut W, n: i32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}

pub(crate) fn write_u32_le<W>(writer: &mut W, n: u32) -> io::Result<()>
where
    W: Write,
{
    let buf = n.to_le_bytes();
    writer.write_all(&buf)
}
