use std::{io, mem};

use bytes::{Buf, BytesMut};

use crate::record::Cigar;

pub fn get_cigar(src: &mut BytesMut, cigar: &mut Cigar, n_cigar_op: usize) -> io::Result<()> {
    let len = mem::size_of::<u32>() * n_cigar_op;

    if src.remaining() < len {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }

    cigar.clear();
    cigar.buf = src.split_to(len);

    Ok(())
}
