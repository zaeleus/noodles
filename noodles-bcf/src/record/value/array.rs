mod values;

pub(crate) use self::values::Values;

pub(crate) enum Array<'a> {
    Int8(Values<'a, i8>),
    Int16(Values<'a, i16>),
    Int32(Values<'a, i32>),
    Float(Values<'a, f32>),
}
