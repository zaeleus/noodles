#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Byte {
    // block_content_id
    External(i32),
    // alphabet, bit_lens
    Huffman(Vec<i32>, Vec<u32>),
}
