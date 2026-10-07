use super::Data;

#[doc(hidden)]
pub enum DataRef<'a> {
    TypedMapEncoded(&'a [u8]),
    Data(Box<dyn Data<'a> + 'a>),
}
