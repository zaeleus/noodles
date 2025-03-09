#![warn(missing_docs)]

//! **noodles-cram** handles the reading and writing of the CRAM format.

#[cfg(feature = "async")]
pub mod r#async;

pub mod codecs;
pub mod container;
pub mod crai;
pub mod file_definition;
pub mod fs;
mod huffman;
pub mod io;
pub mod record;

use md5::{Digest, Md5};

pub use self::{file_definition::FileDefinition, record::Record};

#[deprecated(since = "0.78.0", note = "Use `cram::container` instead.")]
pub use self::container as data_container;

#[deprecated(since = "0.78.0", note = "Use `cram::io::reader::Container` instead.")]
pub use self::io::reader::Container;

#[deprecated(since = "0.78.0", note = "Use `cram::io::reader::Container` instead.")]
pub use self::io::reader::Container as DataContainer;

#[deprecated(since = "0.76.0", note = "Use `cram::fs::index` instead.")]
pub use self::fs::index;

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Reader` instead.")]
pub use self::r#async::io::Reader as AsyncReader;

#[cfg(feature = "async")]
#[deprecated(since = "0.69.0", note = "Use `cram::r#async::io::Writer` instead.")]
pub use self::r#async::io::Writer as AsyncWriter;

const MAGIC_NUMBER: [u8; 4] = *b"CRAM";
const MD5_OUTPUT_SIZE: usize = 16;

// _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.3.2 "Reference MD5 calculation"
fn calculate_normalized_sequence_digest(mut sequence: &[u8]) -> [u8; MD5_OUTPUT_SIZE] {
    const MD5_BLOCK_SIZE: usize = 64;
    const CHUNK_SIZE: usize = 8 * MD5_BLOCK_SIZE;

    let mut hasher = Md5::new();
    let mut buf = [0; CHUNK_SIZE];

    while !sequence.is_empty() {
        let n = sequence
            .iter()
            .position(|b| !b.is_ascii_graphic() || b.is_ascii_lowercase())
            .unwrap_or(sequence.len());

        hasher.update(&sequence[..n]);
        sequence = &sequence[n..];

        // "All lowercase characters are converted to uppercase."
        loop {
            let mut n = 0;

            for (src, dst) in sequence
                .iter()
                .take_while(|b| b.is_ascii_lowercase())
                .zip(&mut buf)
            {
                *dst = src.to_ascii_uppercase();
                n += 1;
            }

            if n == 0 {
                break;
            }

            hasher.update(&buf[..n]);
            sequence = &sequence[n..];
        }

        // "All characters outside of the inclusive range 33 ('!') to 126 ('~') are stripped out."
        let n = sequence
            .iter()
            .position(|b| b.is_ascii_graphic())
            .unwrap_or(sequence.len());

        sequence = &sequence[n..];
    }

    hasher.finalize().into()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_calculate_normalized_sequence_digest() {
        assert_eq!(
            calculate_normalized_sequence_digest(b"ACGT"),
            [
                0xf1, 0xf8, 0xf4, 0xbf, 0x41, 0x3b, 0x16, 0xad, 0x13, 0x57, 0x22, 0xaa, 0x45, 0x91,
                0x04, 0x3e
            ]
        );

        assert_eq!(
            calculate_normalized_sequence_digest(b" AC\tgt\n"),
            [
                0xf1, 0xf8, 0xf4, 0xbf, 0x41, 0x3b, 0x16, 0xad, 0x13, 0x57, 0x22, 0xaa, 0x45, 0x91,
                0x04, 0x3e
            ]
        );

        // _Sequence Alignment/Map Format Specification_ (2021-06-03) § 1.3.2 "Reference MD5
        // calculation"
        assert_eq!(
            calculate_normalized_sequence_digest(b"ACGTACGTACGTACGTACGTACGT...12345!!!"),
            [
                0xdf, 0xab, 0xdb, 0xb3, 0x6e, 0x23, 0x9a, 0x6d, 0xa8, 0x89, 0x57, 0x84, 0x1f, 0x32,
                0xb8, 0xe4
            ]
        );
    }
}
