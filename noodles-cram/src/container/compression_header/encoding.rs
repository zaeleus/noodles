mod kind;

pub use self::kind::Kind;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Encoding {
    kind: Kind,
    args: Vec<u8>,
}

impl Encoding {
    pub fn new(kind: Kind, args: Vec<u8>) -> Self {
        Self { kind, args }
    }

    pub fn kind(&self) -> Kind {
        self.kind
    }

    pub fn args(&self) -> &[u8] {
        &self.args
    }
}
