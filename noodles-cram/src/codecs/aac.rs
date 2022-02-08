#![allow(dead_code)]

mod flags;
mod model;
mod range_coder;

pub use self::{model::Model, range_coder::RangeCoder};

use std::io::{self, Read};

use byteorder::ReadBytesExt;

use self::flags::Flags;

use super::rans_nx16::{decode_pack, decode_pack_meta};

fn arith_decode<R>(reader: &mut R, mut len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let flags = reader.read_u8().map(Flags::from)?;

    if flags.contains(Flags::STRIPE) {
        todo!("arith_decode: decode_stripe");
    }

    let mut p = None;
    let mut n_sym = None;
    let pack_len = len;

    if flags.contains(Flags::PACK) {
        let (q, n, new_len) = decode_pack_meta(reader)?;
        p = Some(q);
        n_sym = Some(n);
        len = new_len;
    }

    let mut data = vec![0; len];

    if flags.contains(Flags::CAT) {
        reader.read_exact(&mut data)?;
    } else if flags.contains(Flags::EXT) {
        todo!("arith_decode: decode_ext");
    } else if flags.contains(Flags::RLE) {
        if flags.contains(Flags::ORDER) {
            todo!("arith_decode: decode_rle_1");
        } else {
            todo!("arith_decode: decode_rle_0");
        }
    } else if flags.contains(Flags::ORDER) {
        todo!("arith_decode: decode_order_1");
    } else {
        todo!("arith_decode: decode_order_0");
    }

    if flags.contains(Flags::PACK) {
        let p = p.unwrap();
        let n_sym = n_sym.unwrap();
        data = decode_pack(&data, &p, n_sym, pack_len)?;
    }

    Ok(data)
}
