mod values;

use std::fmt;

pub use self::values::Values;

pub enum Array<'a> {
    Int8(Box<dyn Values<'a, i8> + 'a>),
    Int16(Box<dyn Values<'a, i16> + 'a>),
    Int32(Box<dyn Values<'a, i32> + 'a>),
    Float(Box<dyn Values<'a, f32> + 'a>),
}

impl fmt::Debug for Array<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Int8(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Int16(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Int32(values) => f.debug_list().entries(values.iter()).finish(),
            Self::Float(values) => f.debug_list().entries(values.iter()).finish(),
        }
    }
}
