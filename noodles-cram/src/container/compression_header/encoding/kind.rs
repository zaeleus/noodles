#[derive(Clone, Copy, Debug, Eq, PartialEq)]
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
