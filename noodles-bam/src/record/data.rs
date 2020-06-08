mod field;
mod reader;
mod value;

pub use self::{field::Field, reader::Reader, value::Value};

use std::ops::Deref;

use self::reader::Fields;

#[derive(Debug)]
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub fn new(bytes: &[u8]) -> Data<'_> {
        Data(bytes)
    }

    pub fn fields(&self) -> Fields<&[u8]> {
        let reader = Reader::new(self.0);
        reader.fields()
    }
}

impl<'a> Deref for Data<'a> {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        self.0
    }
}
