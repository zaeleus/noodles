#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Char(char),
    Int8(i8),
    UInt8(u8),
    Int16(i16),
    UInt16(u16),
    Int32(i32),
    UInt32(u32),
    Float(f32),
    String(String),
    Hex(String),
    Int8Array(Vec<i8>),
    UInt8Array(Vec<u8>),
    Int16Array(Vec<i16>),
    UInt16Array(Vec<u16>),
    Int32Array(Vec<i32>),
    UInt32Array(Vec<u32>),
    FloatArray(Vec<f32>),
}

impl Value {
    pub fn ty(&self) -> char {
        match *self {
            Self::Char(_) => 'A',
            Self::Int8(_) => 'c',
            Self::UInt8(_) => 'C',
            Self::Int16(_) => 's',
            Self::UInt16(_) => 'S',
            Self::Int32(_) => 'i',
            Self::UInt32(_) => 'I',
            Self::Float(_) => 'f',
            Self::String(_) => 'Z',
            Self::Hex(_) => 'H',
            Self::Int8Array(_) => 'B',
            Self::UInt8Array(_) => 'B',
            Self::Int16Array(_) => 'B',
            Self::UInt16Array(_) => 'B',
            Self::Int32Array(_) => 'B',
            Self::UInt32Array(_) => 'B',
            Self::FloatArray(_) => 'B',
        }
    }

    pub fn subtype(&self) -> Option<char> {
        match *self {
            Self::Int8Array(_) => Some('c'),
            Self::UInt8Array(_) => Some('C'),
            Self::Int16Array(_) => Some('s'),
            Self::UInt16Array(_) => Some('S'),
            Self::Int32Array(_) => Some('i'),
            Self::UInt32Array(_) => Some('I'),
            Self::FloatArray(_) => Some('f'),
            _ => None,
        }
    }

    pub fn as_char(&self) -> Option<char> {
        match *self {
            Self::Char(c) => Some(c),
            _ => None,
        }
    }

    pub fn is_char(&self) -> bool {
        self.as_char().is_some()
    }

    pub fn as_int8(&self) -> Option<i8> {
        match *self {
            Self::Int8(b) => Some(b),
            _ => None,
        }
    }

    pub fn is_int8(&self) -> bool {
        self.as_int8().is_some()
    }

    pub fn as_uint8(&self) -> Option<u8> {
        match *self {
            Self::UInt8(b) => Some(b),
            _ => None,
        }
    }

    pub fn is_uint8(&self) -> bool {
        self.as_uint8().is_some()
    }

    pub fn as_int16(&self) -> Option<i16> {
        match *self {
            Self::Int16(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_int16(&self) -> bool {
        self.as_int16().is_some()
    }

    pub fn as_uint16(&self) -> Option<u16> {
        match *self {
            Self::UInt16(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_uint16(&self) -> bool {
        self.as_uint16().is_some()
    }

    pub fn as_int32(&self) -> Option<i32> {
        match *self {
            Self::Int32(i) => Some(i),
            _ => None,
        }
    }

    pub fn is_int32(&self) -> bool {
        self.as_int32().is_some()
    }

    pub fn as_uint32(&self) -> Option<u32> {
        match *self {
            Self::UInt32(i) => Some(i),
            _ => None,
        }
    }

    pub fn is_uint32(&self) -> bool {
        self.as_uint32().is_some()
    }

    pub fn as_float(&self) -> Option<f32> {
        match *self {
            Self::Float(s) => Some(s),
            _ => None,
        }
    }

    pub fn is_float(&self) -> bool {
        self.as_float().is_some()
    }

    pub fn as_str(&self) -> Option<&str> {
        match *self {
            Self::String(ref s) => Some(s),
            _ => None,
        }
    }

    pub fn is_str(&self) -> bool {
        self.as_str().is_some()
    }

    pub fn as_hex(&self) -> Option<&str> {
        match *self {
            Self::Hex(ref h) => Some(h),
            _ => None,
        }
    }

    pub fn is_hex(&self) -> bool {
        self.as_hex().is_some()
    }

    pub fn as_int8_array(&self) -> Option<&[i8]> {
        match *self {
            Self::Int8Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int8_array(&self) -> bool {
        self.as_int8_array().is_some()
    }

    pub fn as_uint8_array(&self) -> Option<&[u8]> {
        match *self {
            Self::UInt8Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint8_array(&self) -> bool {
        self.as_uint8_array().is_some()
    }

    pub fn as_int16_array(&self) -> Option<&[i16]> {
        match *self {
            Self::Int16Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int16_array(&self) -> bool {
        self.as_int16_array().is_some()
    }

    pub fn as_uint16_array(&self) -> Option<&[u16]> {
        match *self {
            Self::UInt16Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint16_array(&self) -> bool {
        self.as_uint16_array().is_some()
    }

    pub fn as_int32_array(&self) -> Option<&[i32]> {
        match *self {
            Self::Int32Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_int32_array(&self) -> bool {
        self.as_int32_array().is_some()
    }

    pub fn as_uint32_array(&self) -> Option<&[u32]> {
        match *self {
            Self::UInt32Array(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_uint32_array(&self) -> bool {
        self.as_uint32_array().is_some()
    }

    pub fn as_float_array(&self) -> Option<&[f32]> {
        match *self {
            Self::FloatArray(ref a) => Some(a),
            _ => None,
        }
    }

    pub fn is_float_array(&self) -> bool {
        self.as_float_array().is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ty() {
        assert_eq!(Value::Char('m').ty(), 'A');
        assert_eq!(Value::Int8(0).ty(), 'c');
        assert_eq!(Value::UInt8(0).ty(), 'C');
        assert_eq!(Value::Int16(0).ty(), 's');
        assert_eq!(Value::UInt16(0).ty(), 'S');
        assert_eq!(Value::Int32(0).ty(), 'i');
        assert_eq!(Value::UInt32(0).ty(), 'I');
        assert_eq!(Value::Float(0.0).ty(), 'f');
        assert_eq!(Value::String(String::from("noodles")).ty(), 'Z');
        assert_eq!(Value::Hex(String::from("cafe")).ty(), 'H');
        assert_eq!(Value::Int8Array(vec![0]).ty(), 'B');
        assert_eq!(Value::UInt8Array(vec![0]).ty(), 'B');
        assert_eq!(Value::Int16Array(vec![0]).ty(), 'B');
        assert_eq!(Value::UInt16Array(vec![0]).ty(), 'B');
        assert_eq!(Value::Int32Array(vec![0]).ty(), 'B');
        assert_eq!(Value::UInt32Array(vec![0]).ty(), 'B');
        assert_eq!(Value::FloatArray(vec![0.0]).ty(), 'B');
    }

    #[test]
    fn test_subtype() {
        assert_eq!(Value::Char('m').subtype(), None);
        assert_eq!(Value::Int8(0).subtype(), None);
        assert_eq!(Value::UInt8(0).subtype(), None);
        assert_eq!(Value::Int16(0).subtype(), None);
        assert_eq!(Value::UInt16(0).subtype(), None);
        assert_eq!(Value::Int32(0).subtype(), None);
        assert_eq!(Value::UInt32(0).subtype(), None);
        assert_eq!(Value::Float(0.0).subtype(), None);
        assert_eq!(Value::String(String::from("noodles")).subtype(), None);
        assert_eq!(Value::Hex(String::from("cafe")).subtype(), None);
        assert_eq!(Value::Int8Array(vec![0]).subtype(), Some('c'));
        assert_eq!(Value::UInt8Array(vec![0]).subtype(), Some('C'));
        assert_eq!(Value::Int16Array(vec![0]).subtype(), Some('s'));
        assert_eq!(Value::UInt16Array(vec![0]).subtype(), Some('S'));
        assert_eq!(Value::Int32Array(vec![0]).subtype(), Some('i'));
        assert_eq!(Value::UInt32Array(vec![0]).subtype(), Some('I'));
        assert_eq!(Value::FloatArray(vec![0.0]).subtype(), Some('f'));
    }
}
