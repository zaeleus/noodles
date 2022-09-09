pub mod aac;
pub mod bzip2;
pub mod fqzcomp;
pub mod gzip;
pub mod lzma;
pub mod name_tokenizer;
pub mod rans;
pub mod rans_nx16;

#[derive(Clone, Debug)]
#[allow(dead_code)]
pub enum Encoder {
    Gzip,
    Bzip2,
    Lzma,
    Rans4x8,
    RansNx16,
}
