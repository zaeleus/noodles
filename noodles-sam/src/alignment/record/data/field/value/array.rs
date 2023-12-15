//! Alignment record data field array value.

mod values;

pub use self::values::Values;

/// An alignment record data field array value.
pub enum Array<'a> {
    /// An 8-bit integer array (`B:c`).
    Int8(Box<dyn Values<'a, i8>>),
    /// An 8-bit unsigned integer array (`B:C`).
    UInt8(Box<dyn Values<'a, u8>>),
    /// A 16-bit integer array (`B:s`).
    Int16(Box<dyn Values<'a, i16>>),
    /// A 16-bit unsigned integer array (`B:S`).
    UInt16(Box<dyn Values<'a, u16>>),
    /// A 32-bit integer array (`B:i`).
    Int32(Box<dyn Values<'a, i32>>),
    /// A 32-bit unsigned integer array (`B:I`).
    UInt32(Box<dyn Values<'a, u32>>),
    /// A single-precision floating-point array (`B:f`).
    Float(Box<dyn Values<'a, f32>>),
}
