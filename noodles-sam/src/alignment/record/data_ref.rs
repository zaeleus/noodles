#![expect(dead_code)]

use super::Data;

enum DataRef<'a> {
    FieldEncoded(&'a [u8]),
    Data(Box<dyn Data + 'a>),
}
