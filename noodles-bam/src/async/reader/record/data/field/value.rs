mod subtype;
mod ty;

pub use self::{subtype::read_subtype, ty::read_type};

use std::convert::TryFrom;

use tokio::io::{self, AsyncBufRead, AsyncBufReadExt, AsyncReadExt};

use crate::record::data::field::{
    value::{Subtype, Type},
    Value,
};

pub async fn read_value<R>(reader: &mut R, ty: Type) -> io::Result<Value>
where
    R: AsyncBufRead + Unpin,
{
    match ty {
        Type::Char => reader.read_u8().await.map(char::from).map(Value::Char),
        Type::Int8 => reader.read_i8().await.map(Value::Int8),
        Type::UInt8 => reader.read_u8().await.map(Value::UInt8),
        Type::Int16 => reader.read_i16_le().await.map(Value::Int16),
        Type::UInt16 => reader.read_u16_le().await.map(Value::UInt16),
        Type::Int32 => reader.read_i32_le().await.map(Value::Int32),
        Type::UInt32 => reader.read_u32_le().await.map(Value::UInt32),
        Type::Float => reader.read_f32_le().await.map(Value::Float),
        Type::String => read_string(reader).await.map(Value::String),
        Type::Hex => read_string(reader).await.map(Value::Hex),
        Type::Array => read_array(reader).await,
    }
}

async fn read_string<R>(reader: &mut R) -> io::Result<String>
where
    R: AsyncBufRead + Unpin,
{
    let mut buf = Vec::new();
    reader.read_until(b'\0', &mut buf).await?;
    buf.pop();
    String::from_utf8(buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

async fn read_array<R>(reader: &mut R) -> io::Result<Value>
where
    R: AsyncBufRead + Unpin,
{
    let subtype = read_subtype(reader).await?;

    let len = reader.read_i32_le().await.and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    match subtype {
        Subtype::Int8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_i8().await?;
                buf.push(value);
            }

            Ok(Value::Int8Array(buf))
        }
        Subtype::UInt8 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_u8().await?;
                buf.push(value);
            }

            Ok(Value::UInt8Array(buf))
        }
        Subtype::Int16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_i16_le().await?;
                buf.push(value);
            }

            Ok(Value::Int16Array(buf))
        }
        Subtype::UInt16 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_u16_le().await?;
                buf.push(value);
            }

            Ok(Value::UInt16Array(buf))
        }
        Subtype::Int32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_i32_le().await?;
                buf.push(value);
            }

            Ok(Value::Int32Array(buf))
        }
        Subtype::UInt32 => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_u32_le().await?;
                buf.push(value);
            }

            Ok(Value::UInt32Array(buf))
        }
        Subtype::Float => {
            let mut buf = Vec::with_capacity(len);

            for _ in 0..len {
                let value = reader.read_f32_le().await?;
                buf.push(value);
            }

            Ok(Value::FloatArray(buf))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[tokio::test]
    async fn test_read_value() -> io::Result<()> {
        async fn t(mut data: &[u8], ty: Type, expected: Value) -> io::Result<()> {
            let actual = read_value(&mut data, ty).await?;
            assert_eq!(actual, expected);
            Ok(())
        }

        t(&[b'n'], Type::Char, Value::Char('n')).await?;
        t(&[0x00], Type::Int8, Value::Int8(0)).await?;
        t(&[0x00], Type::UInt8, Value::UInt8(0)).await?;
        t(&[0x00, 0x00], Type::Int16, Value::Int16(0)).await?;
        t(&[0x00, 0x00], Type::UInt16, Value::UInt16(0)).await?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Int32, Value::Int32(0)).await?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::UInt32, Value::UInt32(0)).await?;
        t(&[0x00, 0x00, 0x00, 0x00], Type::Float, Value::Float(0.0)).await?;
        t(
            &[0x6e, 0x64, 0x6c, 0x73, 0x00],
            Type::String,
            Value::String(String::from("ndls")),
        )
        .await?;
        t(
            &[0x63, 0x61, 0x66, 0x65, 0x00],
            Type::Hex,
            Value::Hex(String::from("cafe")),
        )
        .await?;

        t(
            &[b'c', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int8Array(vec![0]),
        )
        .await?;
        t(
            &[b'C', 0x01, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt8Array(vec![0]),
        )
        .await?;
        t(
            &[b's', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int16Array(vec![0]),
        )
        .await?;
        t(
            &[b'S', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt16Array(vec![0]),
        )
        .await?;
        t(
            &[b'i', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::Int32Array(vec![0]),
        )
        .await?;
        t(
            &[b'I', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::UInt32Array(vec![0]),
        )
        .await?;
        t(
            &[b'f', 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
            Type::Array,
            Value::FloatArray(vec![0.0]),
        )
        .await?;

        Ok(())
    }
}
