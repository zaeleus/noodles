use std::convert::TryFrom;

use crate::num::Itf8;

#[derive(Debug)]
pub enum Kind {
    Null,
    External,
    Golomb,
    Huffman,
    ByteArrayLen,
    ByteArrayStop,
    Beta,
    Subexp,
    GolombRice,
    Gamma,
}

impl TryFrom<Itf8> for Kind {
    type Error = ();

    fn try_from(id: Itf8) -> Result<Self, Self::Error> {
        match id {
            0 => Ok(Self::Null),
            1 => Ok(Self::External),
            2 => Ok(Self::Golomb),
            3 => Ok(Self::Huffman),
            4 => Ok(Self::ByteArrayLen),
            5 => Ok(Self::ByteArrayStop),
            6 => Ok(Self::Beta),
            7 => Ok(Self::Subexp),
            8 => Ok(Self::GolombRice),
            9 => Ok(Self::Gamma),
            _ => Err(()),
        }
    }
}

#[derive(Debug)]
pub struct Encoding {
    kind: Kind,
    args: Vec<u8>,
}

impl Encoding {
    pub fn new(kind: Kind, args: Vec<u8>) -> Self {
        Self { kind, args }
    }
}
