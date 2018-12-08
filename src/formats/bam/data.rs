pub use self::field::Field;
pub use self::reader::Reader;
pub use self::value::Value;

pub mod field;
pub mod reader;
pub mod value;

use std::io;
use std::ops::Deref;

#[derive(Debug)]
pub struct Data(Vec<u8>);

impl Data {
    pub fn new(raw_data: Vec<u8>) -> Data {
        Data(raw_data)
    }

    pub fn parse(&self) -> io::Result<Vec<Field>> {
        let mut reader = Reader::new(&self.0[..]);
        reader.fields().collect()
    }

    pub fn extend(&mut self, field: &Field) {
        self.0.extend_from_slice(field.tag().as_bytes());

        self.0.push(field.value().ty() as u8);

        if let Some(subtype) = field.value().subtype() {
            self.0.push(subtype as u8);
        }

        match field.value() {
            Value::String(s) => {
                self.0.extend_from_slice(s.as_bytes());
                self.0.push(b'\0');
            },
            _ => unimplemented!(),
        }
    }
}

impl Deref for Data {
    type Target = [u8];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}
