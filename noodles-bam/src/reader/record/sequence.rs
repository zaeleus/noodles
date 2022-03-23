use std::io;

use bytes::Buf;
use noodles_sam::{self as sam, record::sequence::Base};

pub(super) fn get_sequence<B>(
    buf: &mut B,
    sequence: &mut sam::record::Sequence,
    l_seq: usize,
) -> io::Result<()>
where
    B: Buf,
{
    let seq_len = (l_seq + 1) / 2;

    if buf.remaining() < seq_len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    let seq = buf.take(seq_len);

    sequence.clear();

    for &b in seq.chunk() {
        sequence.push(decode_base(b >> 4));
        sequence.push(decode_base(b));
    }

    sequence.as_mut().truncate(l_seq);

    buf.advance(seq_len);

    Ok(())
}

fn decode_base(n: u8) -> Base {
    match n & 0x0f {
        0 => Base::Eq,
        1 => Base::A,
        2 => Base::C,
        3 => Base::M,
        4 => Base::G,
        5 => Base::R,
        6 => Base::S,
        7 => Base::V,
        8 => Base::T,
        9 => Base::W,
        10 => Base::Y,
        11 => Base::H,
        12 => Base::K,
        13 => Base::D,
        14 => Base::B,
        15 => Base::N,
        _ => unreachable!(),
    }
}
