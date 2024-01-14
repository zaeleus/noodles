use std::io::{self, Write};

use byteorder::WriteBytesExt;

use super::{Flags, Model, RangeCoder};
use crate::io::writer::num::write_uint7;

pub fn encode(mut flags: Flags, src: &[u8]) -> io::Result<Vec<u8>> {
    use crate::codecs::rans_nx16::encode::pack;

    let mut src = src.to_vec();
    let mut dst = Vec::new();

    dst.write_u8(u8::from(flags))?;

    if !flags.contains(Flags::NO_SIZE) {
        let ulen =
            u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        write_uint7(&mut dst, ulen)?;
    }

    if flags.contains(Flags::STRIPE) {
        let buf = encode_stripe(&src)?;
        dst.extend(buf);
        return Ok(dst);
    }

    let mut pack_header = None;

    if flags.contains(Flags::PACK) {
        match pack::encode(&src) {
            Ok((header, buf)) => {
                pack_header = Some(header);
                src = buf;
            }
            Err(e) if e.kind() == io::ErrorKind::InvalidInput => {
                flags.remove(Flags::PACK);
                dst[0] = u8::from(flags);
            }
            Err(e) => return Err(e),
        }
    }

    if let Some(header) = pack_header {
        dst.write_all(&header)?;
    }

    if flags.contains(Flags::CAT) {
        dst.write_all(&src)?;
    } else if flags.contains(Flags::EXT) {
        encode_ext(&src, &mut dst)?;
    } else if flags.contains(Flags::RLE) {
        if flags.contains(Flags::ORDER) {
            encode_rle_1(&src, &mut dst)?;
        } else {
            encode_rle_0(&src, &mut dst)?;
        }
    } else if flags.contains(Flags::ORDER) {
        encode_order_1(&src, &mut dst)?;
    } else {
        encode_order_0(&src, &mut dst)?;
    }

    Ok(dst)
}

fn encode_stripe(src: &[u8]) -> io::Result<Vec<u8>> {
    const N: usize = 4;

    let mut ulens = Vec::with_capacity(N);
    let mut t = Vec::with_capacity(N);

    for j in 0..N {
        let mut ulen = src.len() / N;

        if src.len() % N > j {
            ulen += 1;
        }

        let chunk = vec![0; ulen];

        ulens.push(ulen);
        t.push(chunk);
    }

    let mut x = 0;
    let mut i = 0;

    while i < src.len() {
        for j in 0..N {
            if x < ulens[j] {
                t[j][x] = src[i + j];
            }
        }

        x += 1;
        i += N;
    }

    let mut chunks = vec![Vec::new(); N];

    for (chunk, s) in chunks.iter_mut().zip(t.iter()) {
        *chunk = encode(Flags::empty(), s)?;
    }

    let mut dst = Vec::new();

    dst.write_u8(N as u8)?;

    for chunk in &chunks {
        let clen = chunk.len() as u32;
        write_uint7(&mut dst, clen)?;
    }

    for chunk in &chunks {
        dst.write_all(chunk)?;
    }

    Ok(dst)
}

fn encode_ext(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    use bzip2::write::BzEncoder;

    let mut encoder = BzEncoder::new(dst, Default::default());
    encoder.write_all(src)?;
    encoder.finish()?;

    Ok(())
}

fn encode_rle_0(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    dst.write_u8(max_sym.overflowing_add(1).0)?;

    let mut model_lit = Model::new(max_sym);
    let mut model_run = vec![Model::new(3); 258];

    let mut range_coder = RangeCoder::default();

    let mut i = 0;

    while i < src.len() {
        let sym = src[i];
        model_lit.encode(dst, &mut range_coder, sym)?;

        let mut run = src[i + 1..].iter().position(|&s| s != sym).unwrap_or(0);
        i += run + 1;

        let mut rctx = usize::from(sym);

        let mut part = run.min(3);
        model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
        rctx = 256;
        run -= part;

        while part == 3 {
            part = run.min(3);
            model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
            rctx = 257;
            run -= part;
        }
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}

fn encode_rle_1(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    dst.write_u8(max_sym.overflowing_add(1).0)?;

    let model_lit_count = usize::from(max_sym) + 1;
    let mut model_lit = vec![Model::new(max_sym); model_lit_count];
    let mut model_run = vec![Model::new(3); 258];

    let mut range_coder = RangeCoder::default();

    let mut i = 0;
    let mut last = 0;

    while i < src.len() {
        let sym = src[i];
        model_lit[last].encode(dst, &mut range_coder, sym)?;

        let mut run = src[i + 1..].iter().position(|&s| s != sym).unwrap_or(0);
        i += run + 1;

        let mut rctx = usize::from(sym);
        last = usize::from(sym);

        let mut part = run.min(3);
        model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
        rctx = 256;
        run -= part;

        while part == 3 {
            part = run.min(3);
            model_run[rctx].encode(dst, &mut range_coder, part as u8)?;
            rctx = 257;
            run -= part;
        }
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}

fn encode_order_0(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    dst.write_u8(max_sym.overflowing_add(1).0)?;

    let mut model = Model::new(max_sym);
    let mut range_coder = RangeCoder::default();

    for &sym in src {
        model.encode(dst, &mut range_coder, sym)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}

fn encode_order_1(src: &[u8], dst: &mut Vec<u8>) -> io::Result<()> {
    let max_sym = src.iter().max().copied().unwrap_or(0);
    dst.write_u8(max_sym.overflowing_add(1).0)?;

    let model_count = usize::from(max_sym) + 1;
    let mut models = vec![Model::new(max_sym); model_count];

    let mut range_coder = RangeCoder::default();

    models[0].encode(dst, &mut range_coder, src[0])?;

    for window in src.windows(2) {
        let sym_0 = usize::from(window[0]);
        let sym_1 = window[1];
        models[sym_0].encode(dst, &mut range_coder, sym_1)?;
    }

    range_coder.range_encode_end(dst)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_ext() -> io::Result<()> {
        use crate::codecs::bzip2;

        let actual = encode(Flags::EXT, b"noodles")?;

        let mut expected = vec![0x04, 0x07];

        let compression_level = ::bzip2::Compression::default();
        let data = bzip2::encode(compression_level, b"noodles")?;
        expected.extend(data);

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_stripe() -> io::Result<()> {
        let actual = encode(Flags::STRIPE, b"noodles")?;

        let expected = [
            0x08, 0x07, 0x04, 0x09, 0x09, 0x09, 0x08, 0x00, 0x02, 0x6f, 0x00, 0xff, 0xa7, 0xab,
            0x62, 0x00, 0x00, 0x02, 0x70, 0x00, 0xff, 0x84, 0x92, 0x1b, 0x00, 0x00, 0x02, 0x74,
            0x00, 0xf7, 0x27, 0xdb, 0x24, 0x00, 0x00, 0x01, 0x65, 0x00, 0xfd, 0x77, 0x20, 0xb0,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_0() -> io::Result<()> {
        let actual = encode(Flags::empty(), b"noodles")?;

        let expected = [
            0x00, 0x07, 0x74, 0x00, 0xf4, 0xe5, 0xb7, 0x4e, 0x50, 0x0f, 0x2e, 0x97, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER, b"noodles")?;

        let expected = [
            0x01, 0x07, 0x74, 0x00, 0xf4, 0xe3, 0x83, 0x41, 0xe2, 0x9a, 0xef, 0x53, 0x50, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_cat() -> io::Result<()> {
        let actual = encode(Flags::CAT, b"noodles")?;
        let expected = [0x20, 0x07, 0x6e, 0x6f, 0x6f, 0x64, 0x6c, 0x65, 0x73];
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_encode_rle_with_order_0() -> io::Result<()> {
        let actual = encode(Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x40, 0x0d, 0x74, 0x00, 0xf3, 0x4b, 0x21, 0x10, 0xa8, 0xe3, 0x84, 0xfe, 0x6b, 0x22,
            0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_rle_with_order_1() -> io::Result<()> {
        let actual = encode(Flags::ORDER | Flags::RLE, b"noooooooodles")?;

        let expected = [
            0x41, 0x0d, 0x74, 0x00, 0xf3, 0x4a, 0x89, 0x79, 0xc1, 0xe8, 0xc3, 0xc5, 0x62, 0x31,
            0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_pack() -> io::Result<()> {
        let actual = encode(Flags::CAT | Flags::PACK, b"noodles")?;

        let expected = [
            0xa0, 0x07, 0x06, 0x64, 0x65, 0x6c, 0x6e, 0x6f, 0x73, 0x04, 0x43, 0x04, 0x12, 0x05,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
