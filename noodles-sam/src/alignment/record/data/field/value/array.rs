/// An alignment record data field array value.
pub enum Array<'a> {
    /// An 8-bit integer array (`B:c`).
    Int8(Box<dyn Iterator<Item = i8> + 'a>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8(Box<dyn Iterator<Item = u8> + 'a>),
    /// A 16-bit integer array (`B:s`).
    Int16(Box<dyn Iterator<Item = i16> + 'a>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16(Box<dyn Iterator<Item = u16> + 'a>),
    /// A 32-bit integer array (`B:i`).
    Int32(Box<dyn Iterator<Item = i32> + 'a>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32(Box<dyn Iterator<Item = u32> + 'a>),
    /// A single-precision floating-point array (`B:f`).
    Float(Box<dyn Iterator<Item = f32> + 'a>),
}
