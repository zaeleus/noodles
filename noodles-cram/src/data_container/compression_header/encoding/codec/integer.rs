#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Integer {
    // block_content_id
    External(i32),
    // offset, m
    Golomb(i32, i32),
    // alphabet, bit_lens
    Huffman(Vec<i32>, Vec<u32>),
    // offset, len
    Beta(i32, u32),
    // offset, k
    Subexp(i32, i32),
    // offset, log2_m
    GolombRice(i32, i32),
    // offset
    Gamma(i32),
}
