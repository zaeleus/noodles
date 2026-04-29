use super::Data;

#[doc(hidden)]
pub enum DataRef<'a> {
    FieldEncoded(&'a [u8]),
    Data(Box<dyn Data + 'a>),
}
