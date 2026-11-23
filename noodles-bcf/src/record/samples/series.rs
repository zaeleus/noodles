//! BCF record samples series.

pub mod value;

use std::{borrow::Cow, io, mem, ops::Range, str};

use memchr::memchr;
use noodles_vcf::{
    self as vcf,
    header::record::value::map::format::{self, Number},
    variant::record::samples::{
        keys::key,
        series::{Value, value::Array},
    },
};

use crate::record::{
    codec::value::{Float, Int8, Int16, Int32},
    value::{Type, array::Values, read_type, read_value},
};

/// A BCF record samples series.
pub struct Series<'r> {
    id: usize,
    ty: Type,
    sample_count: usize,
    src: &'r [u8],
}

impl<'r> Series<'r> {
    /// Returns the name.
    pub fn name<'h>(&self, header: &'h vcf::Header) -> io::Result<&'h str> {
        header
            .string_maps()
            .strings()
            .get_index(self.id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid string map ID"))
    }

    fn len(&self) -> usize {
        self.sample_count
    }

    /// Returns the value at the given index.
    pub fn get(&self, header: &vcf::Header, i: usize) -> Option<Option<io::Result<Value<'r>>>> {
        if i >= self.len() {
            return None;
        }

        let name = match self.name(header) {
            Ok(name) => name,
            Err(e) => return Some(Some(Err(e))),
        };

        if name == key::GENOTYPE {
            return get_genotype_value(self.src, self.ty, header, i);
        }

        let (number, ty) = header
            .formats()
            .get(name)
            .map(|format| (format.number(), format.ty()))
            .expect("missing type definition");

        match (number, ty, self.ty) {
            (Number::Count(0), _, _) => todo!("invalid number for type"),

            (_, _, Type::Int8(0) | Type::Int16(0) | Type::Int32(0) | Type::Float(0)) => Some(Some(
                Err(io::Error::new(io::ErrorKind::InvalidData, "invalid length")),
            )),
            (_, format::Type::Character | format::Type::String, Type::String(0)) => Some(None),

            (Number::Count(1), format::Type::Integer, Type::Int8(1)) => get_i8_value(self.src, i),
            (Number::Count(1), format::Type::Integer, Type::Int16(1)) => get_i16_value(self.src, i),
            (Number::Count(1), format::Type::Integer, Type::Int32(1)) => get_i32_value(self.src, i),
            (Number::Count(1), format::Type::Float, Type::Float(1)) => get_f32_value(self.src, i),
            (
                Number::Count(1),
                format::Type::Integer,
                Type::Int8(_) | Type::Int16(_) | Type::Int32(_),
            )
            | (Number::Count(1), format::Type::Float, Type::Float(_)) => Some(Some(Err(
                io::Error::new(io::ErrorKind::InvalidData, "invalid length"),
            ))),

            (Number::Count(1), format::Type::Character, Type::String(len)) => {
                get_char_value(self.src, len, i)
            }
            (Number::Count(1), format::Type::String, Type::String(len)) => {
                get_string_value(self.src, len, i)
            }

            (_, format::Type::Integer, Type::Int8(len)) => get_i8_array_value(self.src, len, i),
            (_, format::Type::Integer, Type::Int16(len)) => get_i16_array_value(self.src, len, i),
            (_, format::Type::Integer, Type::Int32(len)) => get_i32_array_value(self.src, len, i),
            (_, format::Type::Float, Type::Float(len)) => get_f32_array_value(self.src, len, i),
            (_, format::Type::Character, Type::String(len)) => {
                get_char_array_value(self.src, len, i)
            }
            (_, format::Type::String, Type::String(len)) => {
                get_string_array_value(self.src, len, i)
            }

            _ => Some(Some(Err(type_mismatch_error(self.ty, ty)))),
        }
    }
}

impl vcf::variant::record::samples::Series for Series<'_> {
    fn name<'a, 'h: 'a>(&'a self, header: &'h vcf::Header) -> io::Result<&'a str> {
        header
            .string_maps()
            .strings()
            .get_index(self.id)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "invalid string map ID"))
    }

    fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        i: usize,
    ) -> Option<Option<io::Result<Value<'a>>>> {
        self.get(header, i)
    }

    fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
    ) -> Box<dyn Iterator<Item = io::Result<Option<Value<'a>>>> + 'a> {
        Box::new((0..self.len()).map(|i| {
            self.get(header, i)
                .ok_or_else(|| io::Error::from(io::ErrorKind::InvalidData))?
                .transpose()
        }))
    }
}

pub(super) fn read_series<'a>(src: &mut &'a [u8], sample_count: usize) -> io::Result<Series<'a>> {
    fn size_of(ty: Type) -> usize {
        match ty {
            Type::Int8(n) => mem::size_of::<i8>() * n,
            Type::Int16(n) => mem::size_of::<i16>() * n,
            Type::Int32(n) => mem::size_of::<i32>() * n,
            Type::Float(n) => mem::size_of::<f32>() * n,
            Type::String(n) => mem::size_of::<u8>() * n,
        }
    }

    let id = read_string_map_index(src)?;
    let ty = read_type(src)?.expect("invalid type");

    let len = size_of(ty) * sample_count;
    let (buf, rest) = src.split_at(len);

    *src = rest;

    Ok(Series {
        id,
        ty,
        sample_count,
        src: buf,
    })
}

fn read_string_map_index(src: &mut &[u8]) -> io::Result<usize> {
    match read_value(src)?.and_then(|v| v.as_int()) {
        Some(i) => usize::try_from(i).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)),
        None => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid string map index",
        )),
    }
}

fn range<N>(i: usize, len: usize) -> Range<usize> {
    let size = mem::size_of::<N>();
    let start = size * i * len;
    let end = start + size * len;
    start..end
}

fn get_i8_value(src: &[u8], i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    use crate::record::codec::value::Int8;

    let src = src.get(range::<i8>(i, 1))?;

    match Int8::from(src[0] as i8) {
        Int8::Value(n) => Some(Some(Ok(Value::Integer(i32::from(n))))),
        Int8::Missing => Some(None),
        Int8::EndOfVector | Int8::Reserved(_) => {
            Some(Some(Err(io::Error::from(io::ErrorKind::InvalidData))))
        }
    }
}

fn get_i8_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    let src = src.get(range::<i8>(i, len))?;
    let values = Values::<'_, i8>::new(src);

    let mut iter = values.iter();

    match (iter.next(), iter.next()) {
        (Some(Int8::Missing), None | Some(Int8::EndOfVector)) => Some(None),
        (Some(Int8::EndOfVector), _) => Some(Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty value",
        )))),
        _ => Some(Some(Ok(Value::Array(Array::Integer(Box::new(values)))))),
    }
}

fn get_i16_value(src: &[u8], i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    use crate::record::codec::value::Int16;

    let src = src.get(range::<i16>(i, 1))?;

    // SAFETY: `src` is 2 bytes.
    match Int16::from(i16::from_le_bytes(src.try_into().unwrap())) {
        Int16::Value(n) => Some(Some(Ok(Value::Integer(i32::from(n))))),
        Int16::Missing => Some(None),
        Int16::EndOfVector | Int16::Reserved(_) => {
            Some(Some(Err(io::Error::from(io::ErrorKind::InvalidData))))
        }
    }
}

fn get_i16_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    let src = src.get(range::<i16>(i, len))?;
    let values = Values::<'_, i16>::new(src);

    let mut iter = values.iter();

    match (iter.next(), iter.next()) {
        (Some(Int16::Missing), None | Some(Int16::EndOfVector)) => Some(None),
        (Some(Int16::EndOfVector), _) => Some(Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty value",
        )))),
        _ => Some(Some(Ok(Value::Array(Array::Integer(Box::new(values)))))),
    }
}

fn get_i32_value(src: &[u8], i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    use crate::record::codec::value::Int32;

    let src = src.get(range::<i32>(i, 1))?;

    // SAFETY: `src` is 4 bytes.
    match Int32::from(i32::from_le_bytes(src.try_into().unwrap())) {
        Int32::Value(n) => Some(Some(Ok(Value::Integer(n)))),
        Int32::Missing => Some(None),
        Int32::EndOfVector | Int32::Reserved(_) => {
            Some(Some(Err(io::Error::from(io::ErrorKind::InvalidData))))
        }
    }
}

fn get_i32_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    let src = src.get(range::<i32>(i, len))?;
    let values = Values::<'_, i32>::new(src);

    let mut iter = values.iter();

    match (iter.next(), iter.next()) {
        (Some(Int32::Missing), None | Some(Int32::EndOfVector)) => Some(None),
        (Some(Int32::EndOfVector), _) => Some(Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty value",
        )))),
        _ => Some(Some(Ok(Value::Array(Array::Integer(Box::new(values)))))),
    }
}

fn get_f32_value(src: &[u8], i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    use crate::record::codec::value::Float;

    let src = src.get(range::<f32>(i, 1))?;

    // SAFETY: `src` is 4 bytes.
    match Float::from(f32::from_le_bytes(src.try_into().unwrap())) {
        Float::Value(n) => Some(Some(Ok(Value::Float(n)))),
        Float::Missing => Some(None),
        Float::EndOfVector | Float::Reserved(_) => {
            Some(Some(Err(io::Error::from(io::ErrorKind::InvalidData))))
        }
    }
}

fn get_f32_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    let src = src.get(range::<f32>(i, len))?;
    let values = Values::<'_, f32>::new(src);

    let mut iter = values.iter();

    match (iter.next(), iter.next()) {
        (Some(Float::Missing), None | Some(Float::EndOfVector)) => Some(None),
        (Some(Float::EndOfVector), _) => Some(Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty value",
        )))),
        _ => Some(Some(Ok(Value::Array(Array::Float(Box::new(values)))))),
    }
}

fn get_string(src: &[u8], len: usize, i: usize) -> Option<io::Result<&str>> {
    const NUL: u8 = 0x00;

    let src = src.get(range::<u8>(i, len))?;

    let src = match memchr(NUL, src) {
        Some(i) => &src[..i],
        None => src,
    };

    Some(str::from_utf8(src).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e)))
}

fn get_nonempty_string(src: &[u8], len: usize, i: usize) -> Option<io::Result<&str>> {
    const EMPTY: &str = "";

    match get_string(src, len, i)? {
        Ok(EMPTY) => Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "empty value",
        ))),
        result => Some(result),
    }
}

fn get_char_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    const MISSING: char = '.';

    let s = match get_nonempty_string(src, len, i)? {
        Ok(t) => t,
        Err(e) => return Some(Some(Err(e))),
    };

    let mut chars = s.chars();

    match (chars.next(), chars.next()) {
        (Some(MISSING), None) => Some(None),
        (Some(c), None) => Some(Some(Ok(Value::Character(c)))),
        _ => Some(Some(Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid character",
        )))),
    }
}

fn get_char_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    const MISSING: &str = ".";

    match get_nonempty_string(src, len, i)? {
        Ok(MISSING) => Some(None),
        Ok(s) => Some(Some(Ok(Value::Array(Array::Character(Box::new(s)))))),
        Err(e) => Some(Some(Err(e))),
    }
}

fn get_string_value(src: &[u8], len: usize, i: usize) -> Option<Option<io::Result<Value<'_>>>> {
    const MISSING: &str = ".";

    match get_nonempty_string(src, len, i)? {
        Ok(MISSING) => Some(None),
        Ok(s) => Some(Some(Ok(Value::String(Cow::from(s))))),
        Err(e) => Some(Some(Err(e))),
    }
}

fn get_string_array_value(
    src: &[u8],
    len: usize,
    i: usize,
) -> Option<Option<io::Result<Value<'_>>>> {
    const MISSING: &str = ".";

    match get_nonempty_string(src, len, i)? {
        Ok(MISSING) => Some(None),
        Ok(s) => Some(Some(Ok(Value::Array(Array::String(Box::new(s)))))),
        Err(e) => Some(Some(Err(e))),
    }
}

fn get_genotype_value<'r>(
    src: &'r [u8],
    ty: Type,
    header: &vcf::Header,
    i: usize,
) -> Option<Option<io::Result<Value<'r>>>> {
    use self::value::{Genotype, genotype::Values};
    use crate::record::value::array;

    let values = match ty {
        Type::Int8(len) => src
            .get(range::<i8>(i, len))
            .map(|buf| Values::Int8(array::Values::new(buf)))?,
        Type::Int16(len) => src
            .get(range::<i16>(i, len))
            .map(|buf| Values::Int16(array::Values::new(buf)))?,
        Type::Int32(len) => src
            .get(range::<i32>(i, len))
            .map(|buf| Values::Int32(array::Values::new(buf)))?,
        _ => {
            return Some(Some(Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "invalid type for genotype (GT) value",
            ))));
        }
    };

    Some(Some(Ok(Value::Genotype(Box::new(Genotype::new(
        header.file_format(),
        values,
    ))))))
}

fn type_mismatch_error(actual: Type, expected: format::Type) -> io::Error {
    let resolved_actual = match actual {
        Type::Int8(_) | Type::Int16(_) | Type::Int32(_) => format::Type::Integer,
        Type::Float(_) => format::Type::Float,
        Type::String(_) => format::Type::String,
    };

    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("type mismatch: expected {expected:?}, got {resolved_actual:?}"),
    )
}

#[cfg(test)]
mod tests {
    use noodles_vcf::{
        header::{StringMaps, record::value::Map},
        variant::record::samples::Series as _,
    };

    use super::*;

    const NAME: &str = "value";

    fn build_header_with_format<I>(name: I, number: Number, ty: format::Type) -> vcf::Header
    where
        I: Into<String>,
    {
        let mut header = vcf::Header::builder()
            .add_format(
                name,
                Map::builder()
                    .set_number(number)
                    .set_type(ty)
                    .set_description("")
                    .build()
                    .unwrap(),
            )
            .build();

        *header.string_maps_mut() = StringMaps::try_from(&header).unwrap();

        header
    }

    #[test]
    fn test_get_with_int8_values() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<i32>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Integer(n)) => Some(n),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, // Some(5)
            0x08, // Some(5)
            0x80, // None
            0x81, // EndOfVector
            0x82, // Reserved(0x82)
        ];

        let series = Series {
            id,
            ty: Type::Int8(1),
            sample_count: 5,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        for i in [3, 4] {
            assert!(matches!(
                series.get(&header, i),
                Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        assert!(series.get(&header, 5).is_none());
    }

    #[test]
    fn test_get_with_int8_array_values() {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<i32>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::Integer(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual.as_deref(), expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, 0x08, // Some([Some(5), Some(8)])
            0x0d, 0x80, // Some([Some(13), None])
            0x15, 0x81, // Some([Some(21)])
            0x80, 0x81, // None
        ];

        let series = Series {
            id,
            ty: Type::Int8(2),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some(5), Some(8)]));
        t(&series, &header, 1, Some(&[Some(13), None]));
        t(&series, &header, 2, Some(&[Some(21)]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::Int8(2),
            sample_count: 1,
            src: &[0x81, 0x81],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_int16_values() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<i32>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Integer(n)) => Some(n),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, 0x00, // Some(5)
            0x08, 0x00, // Some(8)
            0x00, 0x80, // None
            0x01, 0x80, // EndOfVector
            0x02, 0x80, // Reserved(0x8002)
        ];

        let series = Series {
            id,
            ty: Type::Int16(1),
            sample_count: 5,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        for i in [3, 4] {
            assert!(matches!(
                series.get(&header, i),
                Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        assert!(series.get(&header, 5).is_none());
    }

    #[test]
    fn test_get_with_int16_array_values() {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<i32>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::Integer(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual.as_deref(), expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, 0x00, 0x08, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x01, 0x80, // Some([Some(21)])
            0x00, 0x80, 0x01, 0x80, // None
        ];

        let series = Series {
            id,
            ty: Type::Int16(2),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some(5), Some(8)]));
        t(&series, &header, 1, Some(&[Some(13), None]));
        t(&series, &header, 2, Some(&[Some(21)]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::Int16(2),
            sample_count: 1,
            src: &[0x01, 0x80, 0x01, 0x80],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_int32_values() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<i32>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Integer(n)) => Some(n),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, 0x00, 0x00, 0x00, // Some(5)
            0x08, 0x00, 0x00, 0x00, // Some(8)
            0x00, 0x00, 0x00, 0x80, // None
            0x01, 0x00, 0x00, 0x80, // EndOfVector
            0x02, 0x00, 0x00, 0x80, // Reserved(0x80000002)
        ];

        let series = Series {
            id,
            ty: Type::Int32(1),
            sample_count: 5,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        for i in [3, 4] {
            assert!(matches!(
                series.get(&header, i),
                Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        assert!(series.get(&header, 5).is_none());
    }

    #[test]
    fn test_get_with_int32_array_values() {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<i32>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::Integer(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual.as_deref(), expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00, // Some([Some(5), Some(8)])
            0x0d, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, // Some([Some(13), None])
            0x15, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // Some([Some(21)])
            0x00, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80, // None
        ];

        let series = Series {
            id,
            ty: Type::Int32(2),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some(5), Some(8)]));
        t(&series, &header, 1, Some(&[Some(13), None]));
        t(&series, &header, 2, Some(&[Some(21)]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::Int32(2),
            sample_count: 1,
            src: &[0x01, 0x00, 0x00, 0x80, 0x01, 0x00, 0x00, 0x80],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_float_values() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<f32>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Float(n)) => Some(n),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Float);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x00, 0x00, 0x00, 0x00, // Some(0.0)
            0x00, 0x00, 0x80, 0x3f, // Some(1.0)
            0x01, 0x00, 0x80, 0x7f, // None
            0x02, 0x00, 0x80, 0x7f, // EndOfVector
            0x03, 0x00, 0x80, 0x7f, // Reserved(0x7f800003)
        ];

        let series = Series {
            id,
            ty: Type::Float(1),
            sample_count: 5,
            src,
        };

        t(&series, &header, 0, Some(0.0));
        t(&series, &header, 1, Some(1.0));
        t(&series, &header, 2, None);

        for i in [3, 4] {
            assert!(matches!(
                series.get(&header, i),
                Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
            ));
        }

        assert!(series.get(&header, 5).is_none());
    }

    #[test]
    fn test_get_with_float_array_values() {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<f32>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::Float(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual.as_deref(), expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Float);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f, // Some([Some(0.0), Some(1.0)])
            0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f, // Some([Some(0.0), None])
            0x00, 0x00, 0x00, 0x00, 0x02, 0x00, 0x80, 0x7f, // Some([Some(0.0)])
            0x01, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f, // None
        ];

        let series = Series {
            id,
            ty: Type::Float(2),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some(0.0), Some(1.0)]));
        t(&series, &header, 1, Some(&[Some(0.0), None]));
        t(&series, &header, 2, Some(&[Some(0.0)]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::Float(2),
            sample_count: 1,
            src: &[0x02, 0x00, 0x80, 0x7f, 0x02, 0x00, 0x80, 0x7f],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_character_value() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<char>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Character(c)) => Some(c),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Character);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();

        #[allow(clippy::byte_char_slices)]
        let src = &[
            b'n', // Some('n')
            b'd', // Some('d')
            b'.', // None
        ];

        let series = Series {
            id,
            ty: Type::String(1),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some('n'));
        t(&series, &header, 1, Some('d'));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());

        let series = Series {
            id,
            ty: Type::String(1),
            sample_count: 1,
            src: &[0x00],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));

        let series = Series {
            id,
            ty: Type::String(2),
            sample_count: 1,
            src: b"nd",
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_character_array_value() {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<char>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::Character(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual.as_deref(), expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Character);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            b'n', b',', b'd', // [Some('n'), Some('d')]
            b'n', b',', b'.', // [Some('n'), None]
            b'n', 0x00, 0x00, // [Some('n')]
            b'.', 0x00, 0x00, // None
        ];

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some('n'), Some('d')]));
        t(&series, &header, 1, Some(&[Some('n'), None]));
        t(&series, &header, 2, Some(&[Some('n')]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 1,
            src: &[0x00, 0x00, 0x00],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 1,
            src: &[b'n', b',', 0x00],
        };
        match series.get(&header, 0).unwrap().transpose().unwrap() {
            Some(Value::Array(Array::Character(values))) => {
                assert!(matches!(
                    values.iter().collect::<Result<Vec<_>, _>>(),
                    Err(e) if e.kind() == io::ErrorKind::InvalidData
                ));
            }
            _ => panic!(),
        }

        let series = Series {
            id,
            ty: Type::String(4),
            sample_count: 1,
            src: b"n,dl",
        };
        match series.get(&header, 0).unwrap().transpose().unwrap() {
            Some(Value::Array(Array::Character(values))) => {
                assert!(matches!(
                    values.iter().collect::<Result<Vec<_>, _>>(),
                    Err(e) if e.kind() == io::ErrorKind::InvalidData
                ));
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_get_with_string_value() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: Option<&str>) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::String(s)) => Some(s),
                Some(_) => panic!(),
                None => None,
            };

            assert_eq!(actual, expected.map(Cow::from));
        }

        let header = build_header_with_format(NAME, Number::Count(1), format::Type::String);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            b'n', 0x00, 0x00, 0x00, // Some("n")
            b'n', b'd', b'l', b's', // Some("ndls")
            b'.', 0x00, 0x00, 0x00, // None
        ];

        let series = Series {
            id,
            ty: Type::String(4),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some("n"));
        t(&series, &header, 1, Some("ndls"));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());

        let series = Series {
            id,
            ty: Type::String(1),
            sample_count: 1,
            src: &[0x00],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));
    }

    #[test]
    fn test_get_with_string_array_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: Option<&[Option<&str>]>,
        ) {
            let actual = match series.get(header, i).unwrap().transpose().unwrap() {
                Some(Value::Array(Array::String(values))) => {
                    Some(values.iter().collect::<Result<Vec<_>, _>>().unwrap())
                }
                Some(_) => panic!(),
                None => None,
            };

            let expected = expected.map(|vs| {
                vs.iter()
                    .map(|value| value.map(Cow::from))
                    .collect::<Vec<_>>()
            });

            assert_eq!(actual, expected);
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::String);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            b'n', 0x00, 0x00, 0x00, // Some([Some("n")])
            b'n', b',', b'l', 0x00, // Some([Some("n"), Some("l")])
            b'n', b',', b'l', b's', // Some([Some("n"), Some("ls")])
            b'.', 0x00, 0x00, 0x00, // None
        ];

        let series = Series {
            id,
            ty: Type::String(4),
            sample_count: 4,
            src,
        };

        t(&series, &header, 0, Some(&[Some("n")]));
        t(&series, &header, 1, Some(&[Some("n"), Some("l")]));
        t(&series, &header, 2, Some(&[Some("n"), Some("ls")]));
        t(&series, &header, 3, None);

        assert!(series.get(&header, 4).is_none());

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 1,
            src: &[0x00, 0x00, 0x00],
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 1,
            src: &[b'n', b',', 0x00],
        };
        match series.get(&header, 0).unwrap().transpose().unwrap() {
            Some(Value::Array(Array::String(values))) => {
                assert!(matches!(
                    values.iter().collect::<Result<Vec<_>, _>>(),
                    Err(e) if e.kind() == io::ErrorKind::InvalidData
                ));
            }
            _ => panic!(),
        }

        Ok(())
    }

    #[test]
    fn test_get_with_genotype_values() -> io::Result<()> {
        use noodles_vcf::variant::record::samples::series::value::genotype::Phasing;

        fn t(
            series: &Series<'_>,
            header: &vcf::Header,
            i: usize,
            expected: &[(Option<usize>, Phasing)],
        ) -> io::Result<()> {
            match series.get(header, i).unwrap().transpose()? {
                Some(Value::Genotype(genotype)) => {
                    let actual = genotype.iter().collect::<io::Result<Vec<_>>>()?;
                    assert_eq!(actual, expected);
                }
                _ => panic!(),
            }

            Ok(())
        }

        let header =
            build_header_with_format(key::GENOTYPE, Number::Count(1), format::Type::String);

        let id = header
            .string_maps()
            .strings()
            .get_index_of(key::GENOTYPE)
            .unwrap();

        // Int8
        let src = &[
            0x02, 0x05, // 0|1
            0x04, 0x81, // 1
        ];

        let series = Series {
            id,
            ty: Type::Int8(2),
            sample_count: 2,
            src,
        };

        t(
            &series,
            &header,
            0,
            &[(Some(0), Phasing::Unphased), (Some(1), Phasing::Phased)],
        )?;
        t(&series, &header, 1, &[(Some(1), Phasing::Unphased)])?;
        assert!(series.get(&header, 2).is_none());

        // Int16
        let src = &[
            0x02, 0x00, 0x01, 0x01, // 0|127
            0x00, 0x01, 0x01, 0x80, // 127
        ];

        let series = Series {
            id,
            ty: Type::Int16(2),
            sample_count: 2,
            src,
        };

        t(
            &series,
            &header,
            0,
            &[(Some(0), Phasing::Unphased), (Some(127), Phasing::Phased)],
        )?;
        t(&series, &header, 1, &[(Some(127), Phasing::Unphased)])?;
        assert!(series.get(&header, 2).is_none());

        // Int32
        let src = &[
            0x02, 0x00, 0x00, 0x00, 0x03, 0x80, 0x00, 0x00, // 0|16384
            0x02, 0x80, 0x00, 0x00, 0x01, 0x00, 0x00, 0x80, // 16384
        ];

        let series = Series {
            id,
            ty: Type::Int32(2),
            sample_count: 2,
            src,
        };

        t(
            &series,
            &header,
            0,
            &[(Some(0), Phasing::Unphased), (Some(16384), Phasing::Phased)],
        )?;
        t(&series, &header, 1, &[(Some(16384), Phasing::Unphased)])?;
        assert!(series.get(&header, 2).is_none());

        // invalid type
        let src = &[0x00, 0x00, 0x00, 0x00];
        let series = Series {
            id,
            ty: Type::Float(1),
            sample_count: 1,
            src,
        };
        assert!(matches!(
            series.get(&header, 0),
            Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_get_with_invalid_scalar_length() {
        for (ty, value_ty) in [
            (format::Type::Integer, Type::Int8(2)),
            (format::Type::Integer, Type::Int16(2)),
            (format::Type::Integer, Type::Int32(2)),
            (format::Type::Float, Type::Float(2)),
        ] {
            let header = build_header_with_format(NAME, Number::Count(1), ty);
            let id = header.string_maps().strings().get_index_of(NAME).unwrap();

            let series = Series {
                id,
                ty: value_ty,
                sample_count: 1,
                src: &[0x00; 8],
            };

            assert!(matches!(
                series.get(&header, 0),
                Some(Some(Err(e))) if e.kind() == io::ErrorKind::InvalidData
            ))
        }
    }

    #[test]
    fn test_get_with_empty_string_type() {
        for (number, ty) in [
            (Number::Count(1), format::Type::Character),
            (Number::Count(2), format::Type::Character),
            (Number::Count(1), format::Type::String),
            (Number::Count(2), format::Type::String),
        ] {
            let header = build_header_with_format(NAME, number, ty);
            let id = header.string_maps().strings().get_index_of(NAME).unwrap();

            let series = Series {
                id,
                ty: Type::String(0),
                sample_count: 1,
                src: &[],
            };

            assert!(matches!(series.get(&header, 0), Some(None)));
            assert!(series.get(&header, 1).is_none());
        }
    }

    #[test]
    fn test_iter() {
        let header = build_header_with_format(NAME, Number::Count(1), format::Type::Integer);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            0x05, // Some(5)
            0x08, // Some(5)
            0x80, // None
        ];

        let series = Series {
            id,
            ty: Type::Int8(1),
            sample_count: 3,
            src,
        };

        let actual: Vec<_> = series
            .iter(&header)
            .map(|result| match result {
                Ok(Some(Value::Integer(n))) => Some(n),
                Ok(None) => None,
                _ => panic!(),
            })
            .collect();

        assert_eq!(actual, [Some(5), Some(8), None]);
    }
}
