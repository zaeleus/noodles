use std::{
    ffi::CString,
    io::{self, Write},
};

use byteorder::{LittleEndian, WriteBytesExt};

use noodles_bam::record::data::field::Value;

pub fn write_value<W>(writer: &mut W, value: &Value) -> io::Result<()>
where
    W: Write,
{
    writer.write_u8(char::from(value.ty()) as u8)?;

    if let Some(subtype) = value.subtype() {
        writer.write_u8(char::from(subtype) as u8)?;
    }

    match value {
        Value::Char(c) => writer.write_u8(*c as u8),
        Value::Int8(n) => writer.write_i8(*n),
        Value::UInt8(n) => writer.write_u8(*n),
        Value::Int16(n) => writer.write_i16::<LittleEndian>(*n),
        Value::UInt16(n) => writer.write_u16::<LittleEndian>(*n),
        Value::Int32(n) => writer.write_i32::<LittleEndian>(*n),
        Value::UInt32(n) => writer.write_u32::<LittleEndian>(*n),
        Value::Float(n) => writer.write_f32::<LittleEndian>(*n),
        Value::String(s) | Value::Hex(s) => {
            let c_str = CString::new(s.as_bytes())
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
            writer.write_all(c_str.as_bytes_with_nul())
        }
        Value::Int8Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_i8(n)?;
            }

            Ok(())
        }
        Value::UInt8Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_u8(n)?;
            }

            Ok(())
        }
        Value::Int16Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_i16::<LittleEndian>(n)?;
            }

            Ok(())
        }
        Value::UInt16Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_u16::<LittleEndian>(n)?;
            }

            Ok(())
        }
        Value::Int32Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_i32::<LittleEndian>(n)?;
            }

            Ok(())
        }
        Value::UInt32Array(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_u32::<LittleEndian>(n)?;
            }

            Ok(())
        }
        Value::FloatArray(values) => {
            writer.write_u32::<LittleEndian>(values.len() as u32)?;

            for &n in values {
                writer.write_f32::<LittleEndian>(n)?;
            }

            Ok(())
        }
    }
}
