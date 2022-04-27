use std::io;

use bytes::{Buf, BytesMut};

use crate::record::Sequence;

pub fn get_sequence(src: &mut BytesMut, sequence: &mut Sequence, l_seq: usize) -> io::Result<()> {
    let seq_len = (l_seq + 1) / 2;

    if src.remaining() < seq_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    sequence.buf = src.split_to(seq_len);
    sequence.len = l_seq;

    Ok(())
}
