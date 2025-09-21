use std::io::{self, Read};

use super::{Flags, Model, RangeCoder};
use crate::io::reader::num::{read_u8, read_uint7_as};

pub fn decode<R>(reader: &mut R, mut len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    use crate::codecs::rans_nx16::decode::pack;

    let flags = read_u8(reader).map(Flags::from)?;

    if !flags.contains(Flags::NO_SIZE) {
        len = read_uint7_as(reader)?;
    }

    if flags.contains(Flags::STRIPE) {
        return decode_stripe(reader, len);
    }

    let mut p = None;
    let mut n_sym = None;
    let pack_len = len;

    if flags.contains(Flags::PACK) {
        let (q, n, new_len) = pack::decode_pack_meta(reader)?;
        p = Some(q);
        n_sym = Some(n);
        len = new_len;
    }

    let mut data = vec![0; len];

    if flags.contains(Flags::CAT) {
        reader.read_exact(&mut data)?;
    } else if flags.contains(Flags::EXT) {
        decode_ext(reader, &mut data)?;
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
        data = pack::decode(&data, &p, n_sym, pack_len)?;
    }

    Ok(data)
}

fn decode_stripe<R>(reader: &mut R, len: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let n = read_u8(reader).map(usize::from)?;
    let mut clens: Vec<usize> = Vec::with_capacity(n);

    for _ in 0..n {
        let clen = read_uint7_as(reader)?;
        clens.push(clen);
    }

    let mut ulens = Vec::with_capacity(n);
    let mut t = Vec::with_capacity(n);

    for j in 0..n {
        let mut ulen = len / n;

        if len % n > j {
            ulen += 1;
        }

        let chunk = decode(reader, ulen)?;

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut dst = vec![0; len];

    for j in 0..n {
        for i in 0..ulens[j] {
            dst[i * n + j] = t[j][i];
        }
    }

    Ok(dst)
}

fn decode_ext<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    use bzip2::read::BzDecoder;

    let mut decoder = BzDecoder::new(reader);
    decoder.read_exact(dst)
}

fn decode_rle_0<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let max_sym = read_u8(reader).map(|n| if n == 0 { u8::MAX } else { n - 1 })?;

    let mut model_lit = Model::new(max_sym);
    let mut model_run = vec![Model::new(3); 258];

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

fn decode_rle_1<R>(reader: &mut R, dst: &mut [u8]) -> io::Result<()>
where
    R: Read,
{
    let max_sym = read_u8(reader).map(|n| if n == 0 { u8::MAX } else { n - 1 })?;

    let mut model_lit = vec![Model::new(max_sym); usize::from(max_sym) + 1];
    let mut model_run = vec![Model::new(3); 258];

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
    let max_sym = read_u8(reader).map(|n| n.overflowing_sub(1).0)?;

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
    let max_sym = read_u8(reader).map(|n| n.overflowing_sub(1).0)?;

    let mut models = vec![Model::new(max_sym); usize::from(max_sym) + 1];

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(reader)?;

    let mut last = 0;

    for b in dst {
        *b = models[last].decode(reader, &mut range_coder)?;
        last = usize::from(*b);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode_order_0() -> io::Result<()> {
        let data = [
            0x00, // flags = {empty}
            0x07, // uncompressed len = 7
            0x74, 0x00, 0xf4, 0xe5, 0xb7, 0x4e, 0x50, 0x0f, 0x2e, 0x97, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_order_1() -> io::Result<()> {
        let data = [
            0x01, // flags = ORDER
            0x07, // uncompressed len = 7
            0x74, 0x00, 0xf4, 0xe3, 0x83, 0x41, 0xe2, 0x9a, 0xef, 0x53, 0x50, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_stripe() -> io::Result<()> {
        let data = [
            0x08, // flags = STRIPE
            0x07, // uncompressed len = 7
            0x04, 0x09, 0x09, 0x09, 0x08, 0x00, 0x02, 0x6f, 0x00, 0xff, 0xa7, 0xab, 0x62, 0x00,
            0x00, 0x02, 0x70, 0x00, 0xff, 0x84, 0x92, 0x1b, 0x00, 0x00, 0x02, 0x74, 0x00, 0xf7,
            0x27, 0xdb, 0x24, 0x00, 0x00, 0x01, 0x65, 0x00, 0xfd, 0x77, 0x20, 0xb0, 0x00, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_uncompressed() -> io::Result<()> {
        let data = [
            0x20, // flags = CAT
            0x07, // uncompressed len = 7
            0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }

    #[test]
    fn test_decode_rle_with_order_0() -> io::Result<()> {
        let data = [
            0x40, // flags = RLE
            0x0d, // uncompressed len = 13
            0x74, 0x00, 0xf3, 0x4b, 0x21, 0x10, 0xa8, 0xe3, 0x84, 0xfe, 0x6b, 0x22, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noooooooodles");

        Ok(())
    }

    #[test]
    fn test_decode_rle_with_order_1() -> io::Result<()> {
        let data = [
            0x41, // flags = ORDER | RLE
            0x0d, // uncompressed len = 13
            0x74, 0x00, 0xf3, 0x4a, 0x89, 0x79, 0xc1, 0xe8, 0xc3, 0xc5, 0x62, 0x31, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noooooooodles");

        Ok(())
    }

    #[test]
    fn test_decode_bit_packing_with_6_symbols() -> io::Result<()> {
        let data = [
            0x80, // flags = PACK
            0x07, // uncompressed len = 7
            0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x44, 0x00, 0xfc, 0x6e, 0x0c, 0xbf,
            0x01, 0xf8, 0x00,
        ];

        let mut reader = &data[..];
        assert_eq!(decode(&mut reader, 0)?, b"noodles");

        Ok(())
    }
}
