use std::{cmp, io};

use super::{
    Models,
    parameters::{Parameters, fqz_decode_params, parameter::Parameter},
};
use crate::{codecs::aac::RangeCoder, io::reader::num::read_uint7_as};

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let uncompressed_size = read_uncompressed_size(&mut src)?;

    let mut params = fqz_decode_params(&mut src)?;

    let mut models = Models::new(params.max_symbol_count, params.selector_count);
    let mut range_coder = RangeCoder::new(&mut src)?;

    let mut i = 0;

    let mut record = Record::default();
    let mut dst = vec![0; uncompressed_size];

    let mut x = 0;
    let mut ctx = 0;

    let mut last_len = 0;
    let mut rev_len = Vec::new();

    while i < uncompressed_size {
        if record.pos == 0 {
            x = fqz_new_record(
                &mut src,
                &mut params,
                &mut range_coder,
                &mut models,
                &mut record,
                last_len,
                &mut rev_len,
            )?;

            last_len = record.len;

            if record.is_duplicate {
                copy_record(&mut dst, i, record.len);

                i += record.len;
                record.pos = 0;

                continue;
            }

            ctx = params.params[x].context;
        }

        let param = &mut params.params[x];
        let q = models.qual[usize::from(ctx)].decode(&mut src, &mut range_coder)?;

        let j = usize::from(q);
        dst[i] = param.quality_map().map(|map| map[j]).unwrap_or(q);

        ctx = fqz_update_context(param, q, &mut record);

        i += 1;
        record.pos -= 1;
    }

    if params.gflags.has_reversed_values() {
        reverse_qualities(&mut dst, uncompressed_size, &rev_len);
    }

    Ok(dst)
}

fn read_uncompressed_size(src: &mut &[u8]) -> io::Result<usize> {
    read_uint7_as(src)
}

#[derive(Debug, Default)]
struct Record {
    rec_no: usize,
    selector: u8,
    len: usize,
    pos: usize,
    is_duplicate: bool,
    q_ctx: u32,
    delta: u32,
    prev_q: u8,
}

fn fqz_new_record(
    src: &mut &[u8],
    parameters: &mut Parameters,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    record: &mut Record,
    mut last_len: usize,
    rev_len: &mut Vec<(bool, usize)>,
) -> io::Result<usize> {
    let mut x = 0;

    if let Some(model) = models.sel.as_mut() {
        record.selector = model.decode(src, range_coder)?;

        if let Some(table) = parameters.selector_table() {
            let i = usize::from(record.selector);
            x = usize::from(table[i]);
        }
    }

    let param = &parameters.params[x];

    if !param.flags().is_fixed_length() || record.rec_no == 0 {
        last_len = read_length(src, range_coder, models)?;
    }

    record.len = last_len;

    if parameters.gflags.has_reversed_values() {
        let rev = models.rev.decode(src, range_coder).map(decode_bool)?;
        let len = record.len;
        rev_len.push((rev, len));
    }

    if param.flags().has_duplicates() {
        record.is_duplicate = models.dup.decode(src, range_coder).map(decode_bool)?;
    }

    record.rec_no += 1;
    record.pos = record.len;
    record.q_ctx = 0;
    record.delta = 0;
    record.prev_q = 0;

    Ok(x)
}

fn fqz_update_context(param: &mut Parameter, q: u8, record: &mut Record) -> u16 {
    let mut ctx = u32::from(param.context);

    let qualities_table = param.qualities_table();
    record.q_ctx = (record.q_ctx << u32::from(param.q_shift))
        .overflowing_add(u32::from(qualities_table[usize::from(q)]))
        .0;

    ctx += (record.q_ctx & ((1 << param.q_bits) - 1)) << param.q_loc;

    if let Some(table) = param.positions_table() {
        let p = cmp::min(record.pos, 1023);
        ctx += u32::from(table[p]) << param.p_loc;
    }

    if let Some(table) = param.deltas_table() {
        let d = cmp::min(record.delta, 255) as usize;
        ctx += u32::from(table[d]) << param.d_loc;

        if record.prev_q != q {
            record.delta += 1;
        }

        record.prev_q = q;
    }

    if param.flags().has_selector() {
        ctx += u32::from(record.selector) << param.s_loc;
    }

    (ctx & 0xffff) as u16
}

fn read_length(
    src: &mut &[u8],
    range_coder: &mut RangeCoder,
    models: &mut Models,
) -> io::Result<usize> {
    let mut buf = [0; 4];

    for (d, model) in buf.iter_mut().zip(&mut models.len) {
        *d = model.decode(src, range_coder)?;
    }

    usize::try_from(u32::from_le_bytes(buf))
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn reverse_qualities(qual: &mut [u8], qual_len: usize, rev_len: &[(bool, usize)]) {
    let mut rec = 0;
    let mut i = 0;

    while i < qual_len {
        let (rev, len) = rev_len[rec];

        if rev {
            let mut j = 0;
            let mut k = len - 1;

            while j < k {
                qual.swap(i + j, i + k);
                j += 1;
                k -= 1;
            }
        }

        i += len;
        rec += 1;
    }
}

fn decode_bool(n: u8) -> bool {
    n != 0
}

fn copy_record(buf: &mut [u8], pos: usize, prev_len: usize) {
    let (src, dst) = buf.split_at_mut(pos);
    let start = src.len() - prev_len;
    let prev_record = &src[start..];
    dst[..prev_len].copy_from_slice(prev_record);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decode() -> io::Result<()> {
        let src = [
            0x19, 0x05, 0x00, 0x00, 0x00, 0x20, 0x03, 0x82, 0x7f, 0x0f, 0x01, 0x01, 0x7d, 0xff,
            0xff, 0x01, 0x84, 0x00, 0x09, 0xff, 0xff, 0xf6, 0x01, 0x65, 0x00, 0x86, 0x2e, 0x98,
            0xea, 0xca, 0x71, 0x6f, 0x22, 0xcd, 0xd8, 0x40,
        ];

        let actual = decode(&src)?;

        let expected = [
            0, 0, 0, 1, 1, 2, 1, 1, 0, 0, // 1
            0, 1, 2, 3, 3, 3, 3, 3, 3, 3, // 2
            2, 1, 1, 0, 0, // 3
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_reverse_qualities() {
        let mut data = Vec::from(b"ndlsndlsndls");
        let len = data.len();
        let rev_len = vec![(false, 4), (true, 4), (false, 4)];

        reverse_qualities(&mut data, len, &rev_len);

        assert_eq!(data, b"ndlssldnndls");
    }

    #[test]
    fn test_decode_bool() {
        assert!(!decode_bool(0));
        assert!(decode_bool(1));
        assert!(decode_bool(2));
        assert!(decode_bool(255));
    }

    #[test]
    fn test_copy_record() {
        let mut buf = *b"ndls\0\0\0\0";
        copy_record(&mut buf[..], 4, 2);
        assert_eq!(buf, *b"ndlsls\0\0");
    }
}
