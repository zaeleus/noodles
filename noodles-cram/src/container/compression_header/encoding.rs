mod kind;

pub use self::kind::Kind;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Encoding {
    Null,
    // block_content_id
    External(i32),
    // offset, m
    Golomb(i32, i32),
    // alphabet, bit_lens
    Huffman(Vec<i32>, Vec<i32>),
    // len_encoding, value_encoding
    ByteArrayLen(Box<Encoding>, Box<Encoding>),
    // stop_byte, block_content_id
    ByteArrayStop(u8, i32),
    // offset, len
    Beta(i32, i32),
    // offset, k
    Subexp(i32, i32),
    // offset, log2_m
    GolombRice(i32, i32),
    // offset
    Gamma(i32),
}
