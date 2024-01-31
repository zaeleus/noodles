use std::io;

use noodles_vcf::{
    header::record::value::map::info::Type,
    variant::record::info::field::{value::Array, Value},
};

use crate::record::{
    codec::value::{Float, Int16, Int32, Int8},
    value::{read_value as read_typed_value, Array as TypedArray},
    Value as TypedValue,
};

pub(super) fn read_value<'a>(src: &mut &'a [u8], ty: Type) -> io::Result<Option<Value<'a>>> {
    match ty {
        Type::Integer => read_integer_value(src),
        Type::Flag => read_flag_value(src),
        Type::Float => read_float_value(src),
        Type::Character => read_character_value(src),
        Type::String => read_string_value(src),
    }
}

fn read_integer_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_typed_value(src)? {
        None
        | Some(TypedValue::Int8(None | Some(Int8::Missing)))
        | Some(TypedValue::Int16(None | Some(Int16::Missing)))
        | Some(TypedValue::Int32(None | Some(Int32::Missing))) => Ok(None),
        Some(TypedValue::Int8(Some(Int8::Value(n)))) => Ok(Some(Value::Integer(i32::from(n)))),
        Some(TypedValue::Int16(Some(Int16::Value(n)))) => Ok(Some(Value::Integer(i32::from(n)))),
        Some(TypedValue::Int32(Some(Int32::Value(n)))) => Ok(Some(Value::Integer(n))),
        Some(TypedValue::Array(TypedArray::Int8(values))) => {
            Ok(Some(Value::Array(Array::Integer(Box::new(values)))))
        }
        Some(TypedValue::Array(TypedArray::Int16(values))) => {
            Ok(Some(Value::Array(Array::Integer(Box::new(values)))))
        }
        Some(TypedValue::Array(TypedArray::Int32(values))) => {
            Ok(Some(Value::Array(Array::Integer(Box::new(values)))))
        }
        v => Err(type_mismatch_error(v, Type::Integer)),
    }
}

fn read_flag_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_typed_value(src)? {
        None | Some(TypedValue::Int8(Some(Int8::Value(1)))) => Ok(Some(Value::Flag)),
        v => Err(type_mismatch_error(v, Type::Flag)),
    }
}

fn read_float_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_typed_value(src)? {
        None | Some(TypedValue::Float(None | Some(Float::Missing))) => Ok(None),
        Some(TypedValue::Float(Some(Float::Value(n)))) => Ok(Some(Value::Float(n))),
        Some(TypedValue::Array(TypedArray::Float(values))) => {
            Ok(Some(Value::Array(Array::Float(Box::new(values)))))
        }
        v => Err(type_mismatch_error(v, Type::Float)),
    }
}

fn read_character_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_typed_value(src)? {
        None | Some(TypedValue::String(None)) => Ok(None),
        Some(TypedValue::String(Some(s))) => match s.len() {
            0 => unreachable!(),
            1 => Ok(Some(Value::Character(s.chars().next().unwrap()))),
            _ => Ok(Some(Value::Array(Array::Character(Box::new(s))))),
        },
        v => Err(type_mismatch_error(v, Type::Character)),
    }
}

fn read_string_value<'a>(src: &mut &'a [u8]) -> io::Result<Option<Value<'a>>> {
    match read_typed_value(src)? {
        None | Some(TypedValue::String(None)) => Ok(None),
        Some(TypedValue::String(Some(s))) => Ok(Some(Value::String(s))),
        v => Err(type_mismatch_error(v, Type::String)),
    }
}

fn type_mismatch_error(value: Option<TypedValue>, expected: Type) -> io::Error {
    let actual = value.map(|v| match v {
        TypedValue::Int8(_) | TypedValue::Int16(_) | TypedValue::Int32(_) => Type::Integer,
        TypedValue::Float(_) => Type::Float,
        TypedValue::String(_) => Type::String,
        TypedValue::Array(TypedArray::Int8(_) | TypedArray::Int16(_) | TypedArray::Int32(_)) => {
            Type::Integer
        }
        TypedValue::Array(TypedArray::Float(_)) => Type::Float,
    });

    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("type mismatch: expected {expected:?}, got {actual:?}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_value_with_integer_value() {
        fn t(mut src: &[u8], expected: Option<i32>) {
            let result = read_value(&mut src, Type::Integer);

            if let Some(n) = expected {
                assert!(matches!(result, Ok(Some(Value::Integer(m))) if m == n));
            } else {
                assert!(matches!(result, Ok(None)));
            }
        }

        // None
        t(&[0x00], None);

        // Some(Value::Int8(None))
        t(&[0x01], None);
        // Some(Value::Int8(Int8::Missing))
        t(&[0x11, 0x80], None);
        // Some(Value::Int8(Some(Int8::Value(8))))
        t(&[0x11, 0x08], Some(8));

        // Some(Value::Int16(None))
        t(&[0x02], None);
        // Some(Value::Int16(Int16::Missing))
        t(&[0x12, 0x00, 0x80], None);
        // Some(Value::Int16(Some(Int16::Value(13))))
        t(&[0x12, 0x0d, 0x00], Some(13));

        // Some(Value::Int32(None))
        t(&[0x03], None);
        // Some(Value::Int32(Int32::Missing))
        t(&[0x13, 0x00, 0x00, 0x00, 0x80], None);
        // Some(Value::Int32(Some(Int32::Value(21))))
        t(&[0x13, 0x15, 0x00, 0x00, 0x00], Some(21));
    }

    #[test]
    fn test_read_value_with_integer_array_value() {
        fn t(mut src: &[u8], expected: &[Option<i32>]) {
            match read_value(&mut src, Type::Integer) {
                Ok(Some(Value::Array(Array::Integer(values)))) => {
                    assert!(matches!(
                        values.iter().collect::<io::Result<Vec<_>>>(),
                        Ok(vs) if vs == expected
                    ));
                }
                _ => panic!(),
            }
        }

        // Some(Value::Array(Array::Int8([Some(8), Some(13)])))
        t(&[0x21, 0x08, 0x0d], &[Some(8), Some(13)]);
        // Some(Value::Array(Array::Int8([Some(8), None])))
        t(&[0x21, 0x08, 0x80], &[Some(8), None]);

        // Some(Value::Array(Array::Int16([Some(21), Some(34)])))
        t(&[0x22, 0x15, 0x00, 0x22, 0x00], &[Some(21), Some(34)]);
        // Some(Value::Array(Array::Int16([Some(21), None])))
        t(&[0x22, 0x15, 0x00, 0x00, 0x80], &[Some(21), None]);

        // Some(Value::Array(Array::Int32([Some(55), Some(89)])))
        t(
            &[0x23, 0x37, 0x00, 0x00, 0x00, 0x59, 0x00, 0x00, 0x00],
            &[Some(55), Some(89)],
        );
        // Some(Value::Array(Array::Int32([Some(55), None])))
        t(
            &[0x23, 0x37, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80],
            &[Some(55), None],
        );
    }

    #[test]
    fn test_read_value_with_flag_value() {
        fn t(mut src: &[u8]) {
            assert!(matches!(
                read_value(&mut src, Type::Flag),
                Ok(Some(Value::Flag))
            ));
        }

        // None
        t(&[0x00]);
        // Some(Value::Int8(Some(Int8::Value(1))))
        t(&[0x11, 0x01]);
    }

    #[test]
    fn test_read_value_with_float_value() {
        fn t(mut src: &[u8], expected: Option<f32>) {
            let result = read_value(&mut src, Type::Float);

            if let Some(n) = expected {
                assert!(matches!(result, Ok(Some(Value::Float(m))) if m == n));
            } else {
                assert!(matches!(result, Ok(None)));
            }
        }

        // None
        t(&[0x00], None);
        // Some(Value::Float(None))
        t(&[0x05], None);
        // Some(Value::Float(Some(Float::Missing)))
        t(&[0x15, 0x01, 0x00, 0x80, 0x7f], None);

        // Some(Value::Float(Some(Float::Value(0.0))))
        t(&[0x15, 0x00, 0x00, 0x00, 0x00], Some(0.0));
    }

    #[test]
    fn test_read_value_with_float_array_value() {
        fn t(mut src: &[u8], expected: &[Option<f32>]) {
            match read_value(&mut src, Type::Float) {
                Ok(Some(Value::Array(Array::Float(values)))) => {
                    assert!(matches!(
                        values.iter().collect::<io::Result<Vec<_>>>(),
                        Ok(vs) if vs == expected
                    ));
                }
                _ => panic!(),
            }
        }

        // Some(Value::Array(Array::Float([0.0, 1.0])))
        t(
            &[0x25, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80, 0x3f],
            &[Some(0.0), Some(1.0)],
        );
        // Some(Value::Array(Array::Float([0.0, None])))
        t(
            &[0x25, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x80, 0x7f],
            &[Some(0.0), None],
        );
    }

    #[test]
    fn test_read_value_with_character_value() {
        fn t(mut src: &[u8], expected: Option<char>) {
            let result = read_value(&mut src, Type::Character);

            if let Some(c) = expected {
                assert!(matches!(result, Ok(Some(Value::Character(d))) if d == c));
            } else {
                assert!(matches!(result, Ok(None)));
            }
        }

        // None
        t(&[0x00], None);
        // Some(Value::String(None))
        t(&[0x07], None);

        // Some(Value::String(Some(String::from("n"))))
        t(&[0x17, 0x6e], Some('n'));
    }

    #[test]
    fn test_read_value_with_character_array_value() {
        fn t(mut src: &[u8], expected: &[Option<char>]) {
            match read_value(&mut src, Type::Character) {
                Ok(Some(Value::Array(Array::Character(values)))) => {
                    assert!(matches!(
                        values.iter().collect::<io::Result<Vec<_>>>(),
                        Ok(vs) if vs == expected
                    ));
                }
                _ => panic!(),
            }
        }

        // Some(Value::String(Some(String::from("n,d"))))
        t(&[0x37, 0x6e, 0x2c, 0x64], &[Some('n'), Some('d')]);
        // Some(Value::String(Some(String::from("n,."))))
        t(&[0x37, 0x6e, 0x2c, 0x2e], &[Some('n'), None]);
    }

    #[test]
    fn test_read_value_with_string_value() {
        fn t(mut src: &[u8], expected: Option<&str>) {
            let result = read_value(&mut src, Type::String);

            if let Some(s) = expected {
                assert!(matches!(result, Ok(Some(Value::String(t))) if t == s));
            } else {
                assert!(matches!(result, Ok(None)));
            }
        }

        // None
        t(&[0x00], None);

        // Some(Value::String(None))
        t(&[0x07], None);
        // Some(Value::String(Some(String::from("ndls"))))
        t(&[0x47, 0x6e, 0x64, 0x6c, 0x73], Some("ndls"));
    }
}
