mod itf8;
mod ltf8;

pub use self::{
    itf8::{read_itf8, write_itf8},
    ltf8::read_ltf8,
};

pub type Itf8 = i32;
pub type Ltf8 = i64;
