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
    Gzip(flate2::Compression),
    Bzip2(::bzip2::Compression),
    Lzma(u32),
    Rans4x8(rans::Order),
    RansNx16(rans_nx16::Flags),
    AdaptiveArithmeticCoding(aac::Flags),
}
