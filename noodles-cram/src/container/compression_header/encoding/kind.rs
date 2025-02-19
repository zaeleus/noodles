#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Kind {
    Null,
    External,
    Golomb,
    Huffman,
    ByteArrayLength,
    ByteArrayStop,
    Beta,
    Subexp,
    GolombRice,
    Gamma,
}
