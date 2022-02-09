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
            decode_rle_1(reader, &mut data)?;
        } else {
            decode_rle_0(reader, &mut data)?;
        }
    } else if flags.contains(Flags::ORDER) {
        decode_order_1(reader, &mut data)?;
    } else {
        decode_order_0(reader, &mut data)?;
    }

    if flags.contains(Flags::PACK) {
        let p = p.unwrap();
        let n_sym = n_sym.unwrap();
        data = decode_pack(&data, &p, n_sym, pack_len)?;
    }

    Ok(data)
}

fn decode_rle_0<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<()>
where
    R: Read,
{
    let max_sym = reader.read_u8().map(|n| if n == 0 { u8::MAX } else { n })?;
    let mut model_lit = Model::new(max_sym);
    let mut model_run = vec![Model::new(4); 258];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(reader)?;

    let mut i = 0;

    while i < dst.len() {
        let b = model_lit.decode(reader, &mut range_coder)?;
        dst[i] = b;

        let mut part = model_run[usize::from(b)].decode(reader, &mut range_coder)?;
        let mut run = usize::from(part);
        let mut rctx = 256;

        while part == 3 {
            part = model_run[rctx].decode(reader, &mut range_coder)?;
            rctx = 257;
            run += usize::from(part);
        }

        for j in 1..=run {
            dst[i + j] = b;
        }

        i += run + 1;
    }

    Ok(())
}

fn decode_rle_1<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<()>
where
    R: Read,
{
    let max_sym = reader.read_u8().map(|n| if n == 0 { u8::MAX } else { n })?;
    let mut model_lit = vec![Model::new(max_sym); usize::from(max_sym)];
    let mut model_run = vec![Model::new(4); 258];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(reader)?;

    let mut i = 0;
    let mut last = 0;

    while i < dst.len() {
        let b = model_lit[last].decode(reader, &mut range_coder)?;
        dst[i] = b;
        last = usize::from(b);

        let mut part = model_run[last].decode(reader, &mut range_coder)?;
        let mut run = usize::from(part);
        let mut rctx = 256;

        while part == 3 {
            part = model_run[rctx].decode(reader, &mut range_coder)?;
            rctx = 257;
            run += usize::from(part);
        }

        for j in 1..=run {
            dst[i + j] = b;
        }

        i += run + 1;
    }

    Ok(())
}

fn decode_order_0<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<()>
where
    R: Read,
{
    let max_sym = reader.read_u8().map(|n| if n == 0 { u8::MAX } else { n })?;
    let mut model = Model::new(max_sym);

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(reader)?;

    for b in dst {
        *b = model.decode(reader, &mut range_coder)?;
    }

    Ok(())
}

fn decode_order_1<R>(reader: &mut R, dst: &mut Vec<u8>) -> io::Result<()>
where
    R: Read,
{
    let max_sym = reader.read_u8().map(|n| if n == 0 { u8::MAX } else { n })?;
    let mut models = vec![Model::new(max_sym); usize::from(max_sym)];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(reader)?;

    let mut last = 0;

    for b in dst {
        *b = models[last].decode(reader, &mut range_coder)?;
        last = usize::from(*b);
    }

    Ok(())
}
