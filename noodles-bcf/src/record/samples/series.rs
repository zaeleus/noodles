//! BCF record samples series.

pub mod value;

use std::{borrow::Cow, io, mem, ops::Range, str};

use noodles_vcf::{
    self as vcf,
    header::record::value::map::format::{self, Number},
    variant::record::samples::{
        keys::key,
        series::{Value, value::Array},
    },
};

use crate::record::value::{Type, array::Values, read_type, read_value};

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
        let name = match self.name(header) {
            Ok(name) => name,
            Err(e) => return Some(Some(Err(e))),
        };

        if name == key::GENOTYPE {
            match self.ty {
                Type::Int8(len) => return get_genotype_value(self.src, header, len, i),
                _ => todo!("unhandled type"),
            }
        }

        let (number, ty) = header
            .formats()
            .get(name)
            .map(|format| (format.number(), format.ty()))
            .expect("missing type definition");

        let value = match (number, ty, self.ty) {
            (Number::Count(0), _, _) => todo!("invalid number for type"),

            (_, _, Type::Int8(0) | Type::Int16(0) | Type::Int32(0) | Type::Float(0)) => {
                return Some(Some(Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid length",
                ))));
            }

            (Number::Count(1), format::Type::Integer, Type::Int8(len)) => {
                get_i8_value(self.src, len, i)
            }
            (Number::Count(1), format::Type::Integer, Type::Int16(len)) => {
                get_i16_value(self.src, len, i)
            }
            (Number::Count(1), format::Type::Integer, Type::Int32(len)) => {
                get_i32_value(self.src, len, i)
            }
            (Number::Count(1), format::Type::Float, Type::Float(len)) => {
                get_f32_value(self.src, len, i)
            }
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

            _ => todo!("unhandled type"),
        };

        match value {
            Some(Some(value)) => Some(Some(Ok(value))),
            Some(None) => Some(None),
            None => None,
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

fn get_i8_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int8;

    let src = src.get(range::<i8>(i, len))?;

    let value = match Int8::from(src[0] as i8) {
        Int8::Value(n) => Some(Value::Integer(i32::from(n))),
        Int8::Missing => None,
        Int8::EndOfVector | Int8::Reserved(_) => todo!(),
    };

    Some(value)
}

fn get_i8_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let src = src.get(range::<i8>(i, len))?;
    let values = Values::<'_, i8>::new(src);
    Some(Some(Value::Array(Array::Integer(Box::new(values)))))
}

fn get_i16_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int16;

    let src = src.get(range::<i16>(i, len))?;

    // SAFETY: `src` is 2 bytes.
    let value = match Int16::from(i16::from_le_bytes(src.try_into().unwrap())) {
        Int16::Value(n) => Some(Value::Integer(i32::from(n))),
        Int16::Missing => None,
        Int16::EndOfVector | Int16::Reserved(_) => todo!(),
    };

    Some(value)
}

fn get_i16_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let src = src.get(range::<i16>(i, len))?;
    let values = Values::<'_, i16>::new(src);
    Some(Some(Value::Array(Array::Integer(Box::new(values)))))
}

fn get_i32_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Int32;

    let src = src.get(range::<i32>(i, len))?;

    // SAFETY: `src` is 2 bytes.
    let value = match Int32::from(i32::from_le_bytes(src.try_into().unwrap())) {
        Int32::Value(n) => Some(Value::Integer(n)),
        Int32::Missing => None,
        Int32::EndOfVector | Int32::Reserved(_) => todo!(),
    };

    Some(value)
}

fn get_i32_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let src = src.get(range::<i32>(i, len))?;
    let values = Values::<'_, i32>::new(src);
    Some(Some(Value::Array(Array::Integer(Box::new(values)))))
}

fn get_f32_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    use crate::record::codec::value::Float;

    let src = src.get(range::<f32>(i, len))?;

    // SAFETY: `src` is 2 bytes.
    let value = match Float::from(f32::from_le_bytes(src.try_into().unwrap())) {
        Float::Value(n) => Some(Value::Float(n)),
        Float::Missing => None,
        Float::EndOfVector | Float::Reserved(_) => todo!(),
    };

    Some(value)
}

fn get_f32_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let src = src.get(range::<f32>(i, len))?;
    let values = Values::<'_, f32>::new(src);
    Some(Some(Value::Array(Array::Float(Box::new(values)))))
}

fn get_string(src: &[u8], len: usize, i: usize) -> Option<&str> {
    const NUL: u8 = 0x00;

    let src = src.get(range::<u8>(i, len))?;

    let src = match src.iter().position(|&b| b == NUL) {
        Some(i) => &src[..i],
        None => src,
    };

    Some(
        str::from_utf8(src)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .unwrap(), // TODO
    )
}

fn get_char_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    const MISSING: char = '.';

    let s = get_string(src, len, i)?;

    // TODO
    let c = s.chars().next().unwrap();

    match c {
        MISSING => Some(None),
        _ => Some(Some(Value::Character(c))),
    }
}

fn get_char_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let s = get_string(src, len, i)?;
    Some(Some(Value::Array(Array::Character(Box::new(s)))))
}

fn get_string_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    const MISSING: &str = ".";

    match get_string(src, len, i)? {
        MISSING => Some(None),
        s => Some(Some(Value::String(Cow::from(s)))),
    }
}

fn get_string_array_value(src: &[u8], len: usize, i: usize) -> Option<Option<Value<'_>>> {
    let s = get_string(src, len, i)?;
    Some(Some(Value::Array(Array::String(Box::new(s)))))
}

fn get_genotype_value<'r>(
    src: &'r [u8],
    header: &vcf::Header,
    len: usize,
    i: usize,
) -> Option<Option<io::Result<Value<'r>>>> {
    use self::value::Genotype;

    let src = src.get(range::<i8>(i, len))?;

    Some(Some(Ok(Value::Genotype(Box::new(Genotype::new(
        header.file_format(),
        src,
    ))))))
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
        ];

        let series = Series {
            id,
            ty: Type::Int8(1),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());
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
        t(&series, &header, 3, Some(&[None])); // FIXME

        assert!(series.get(&header, 4).is_none());
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
        ];

        let series = Series {
            id,
            ty: Type::Int16(1),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());
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
        t(&series, &header, 3, Some(&[None])); // FIXME

        assert!(series.get(&header, 4).is_none());
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
        ];

        let series = Series {
            id,
            ty: Type::Int32(1),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some(5));
        t(&series, &header, 1, Some(8));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());
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
        t(&series, &header, 3, Some(&[None])); // FIXME

        assert!(series.get(&header, 4).is_none());
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
        ];

        let series = Series {
            id,
            ty: Type::Float(1),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, Some(0.0));
        t(&series, &header, 1, Some(1.0));
        t(&series, &header, 2, None);

        assert!(series.get(&header, 3).is_none());
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
        t(&series, &header, 3, Some(&[None])); // FIXME

        assert!(series.get(&header, 4).is_none());
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
    }

    #[test]
    fn test_get_with_character_array_value() {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: &[Option<char>]) {
            match series.get(header, i).unwrap().unwrap().unwrap() {
                Value::Array(Array::Character(values)) => {
                    assert_eq!(
                        values.iter().collect::<Result<Vec<_>, _>>().unwrap(),
                        expected,
                    );
                }
                _ => panic!(),
            }
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::Character);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            b'n', b',', b'd', // [Some('n'), Some('d')]
            b'n', b',', b'.', // [Some('n'), None]
            b'n', 0x00, 0x00, // [Some('n')]
        ];

        let series = Series {
            id,
            ty: Type::String(3),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, &[Some('n'), Some('d')]);
        t(&series, &header, 1, &[Some('n'), None]);
        t(&series, &header, 2, &[Some('n')]);

        assert!(series.get(&header, 3).is_none());
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
    }

    #[test]
    fn test_get_with_string_array_value() -> Result<(), Box<dyn std::error::Error>> {
        fn t(series: &Series<'_>, header: &vcf::Header, i: usize, expected: &[Option<&str>]) {
            let expected: Vec<_> = expected.iter().map(|value| value.map(Cow::from)).collect();

            match series.get(header, i).unwrap().unwrap().unwrap() {
                Value::Array(Array::String(values)) => {
                    assert_eq!(
                        values.iter().collect::<Result<Vec<_>, _>>().unwrap(),
                        expected,
                    );
                }
                _ => panic!(),
            }
        }

        let header = build_header_with_format(NAME, Number::Count(2), format::Type::String);
        let id = header.string_maps().strings().get_index_of(NAME).unwrap();
        let src = &[
            b'n', 0x00, 0x00, 0x00, // "n"
            b'n', b',', b'l', 0x00, // "n,l"
            b'n', b',', b'l', b's', // "n,ls"
        ];

        let series = Series {
            id,
            ty: Type::String(4),
            sample_count: 3,
            src,
        };

        t(&series, &header, 0, &[Some("n")]);
        t(&series, &header, 1, &[Some("n"), Some("l")]);
        t(&series, &header, 2, &[Some("n"), Some("ls")]);

        assert!(series.get(&header, 3).is_none());

        Ok(())
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
