//! CRAM block content codecs.

pub(crate) mod aac;
pub(crate) mod bzip2;
pub(crate) mod fqzcomp;
pub(crate) mod gzip;
pub(crate) mod lzma;
pub(crate) mod name_tokenizer;
pub mod rans_4x8;
pub mod rans_nx16;

/// A CRAM block content encoder.
#[derive(Clone, Debug)]
pub enum Encoder {
    /// gzip
    Gzip(flate2::Compression),
    /// bzip2
    Bzip2(::bzip2::Compression),
    /// xz
    Lzma(u32),
    /// rANS 4x8
    Rans4x8(rans_4x8::Order),
    /// rANS Nx16
    RansNx16(rans_nx16::Flags),
    /// adaptive arithmetic coder
    AdaptiveArithmeticCoding(aac::Flags),
    /// name tokenizer
    NameTokenizer,
    /// fqzcomp
    Fqzcomp,
}
