#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Type {
    Int8(usize),
    Int16(usize),
    Int32(usize),
    Float(usize),
    String(usize),
}
