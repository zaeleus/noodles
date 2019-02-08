pub use self::field::Field;
pub use self::reader::Reader;
pub use self::value::Value;

pub mod field;
pub mod reader;
pub mod value;

use std::ops::Deref;

use self::reader::Fields;

#[derive(Debug)]
pub struct Data<'a>(&'a [u8]);

impl<'a> Data<'a> {
    pub fn new(bytes: &[u8]) -> Data {
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
