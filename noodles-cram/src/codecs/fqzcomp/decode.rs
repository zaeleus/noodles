use std::{cmp, io};

use super::{
    Models,
    parameter::{self, Parameter},
    parameters::{self, Parameters, fqz_decode_params},
};
use crate::{codecs::aac::RangeCoder, io::reader::num::read_uint7_as};

pub fn decode(mut src: &[u8]) -> io::Result<Vec<u8>> {
    let buf_len = read_uint7_as(&mut src)?;

    let mut params = fqz_decode_params(&mut src)?;

    let mut range_coder = RangeCoder::default();
    range_coder.range_decode_create(&mut src)?;

    let mut models = Models::new(params.max_sym, params.max_sel);

    let mut i = 0;

    let mut record = Record::default();
    let mut dst = vec![0; buf_len];

    let mut x = 0;
    let mut ctx = 0;

    let mut last_len = 0;
    let mut rev_len = Vec::new();

    while i < buf_len {
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

            last_len = record.rec_len;

            if record.is_dup {
                for j in 0..record.rec_len {
                    dst[i + j] = dst[i + j - record.rec_len];
                }

                i += record.rec_len;
                record.pos = 0;

                continue;
            }

            ctx = params.params[x].context;
        }

        let param = &mut params.params[x];
        let q = models.qual[usize::from(ctx)].decode(&mut src, &mut range_coder)?;

        dst[i] = if let Some(q_map) = param.q_map.as_deref() {
            q_map[usize::from(q)]
        } else {
            q
        };

        ctx = fqz_update_context(param, q, &mut record);

        i += 1;
        record.pos -= 1;
    }

    if params.gflags.contains(parameters::Flags::DO_REV) {
        reverse_qualities(&mut dst, buf_len, &rev_len);
    }

    Ok(dst)
}

#[derive(Debug, Default)]
struct Record {
    rec: usize,
    sel: u8,
    rec_len: usize,
    pos: usize,
    is_dup: bool,
    qctx: u32,
    delta: u32,
    prevq: u8,
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
    let mut sel = 0;
    let mut x = 0;

    if parameters.max_sel > 0 {
        sel = models.sel.decode(src, range_coder)?;

        if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
            x = usize::from(parameters.s_tab[usize::from(sel)]);
        }
    }

    record.sel = sel;

    let param = &parameters.params[x];

    let is_fixed_len = param.flags.contains(parameter::Flags::DO_LEN);

    if !is_fixed_len || record.rec == 0 {
        last_len = decode_length(src, range_coder, models)? as usize;
    }

    record.rec_len = last_len;
    record.pos = record.rec_len;

    if parameters.gflags.contains(parameters::Flags::DO_REV) {
        let rev = models.rev.decode(src, range_coder).map(|n| n == 1)?;
        let len = record.rec_len;
        rev_len.push((rev, len));
    }

    record.rec += 1;

    if param.flags.contains(parameter::Flags::DO_DEDUP) {
        record.is_dup = models.dup.decode(src, range_coder)? == 1;
    }

    record.qctx = 0;
    record.delta = 0;
    record.prevq = 0;

    Ok(x)
}

fn fqz_update_context(param: &mut Parameter, q: u8, record: &mut Record) -> u16 {
    use parameter::Flags;

    let mut ctx = u32::from(param.context);

    record.qctx = (record.qctx << u32::from(param.q_shift))
        .overflowing_add(u32::from(param.q_tab[usize::from(q)]))
        .0;

    ctx += (record.qctx & ((1 << param.q_bits) - 1)) << param.q_loc;

    if let Some(p_tab) = param.p_tab.as_deref() {
        let p = cmp::min(record.pos, 1023);
        ctx += u32::from(p_tab[p]) << param.p_loc;
    }

    if let Some(d_tab) = param.d_tab.as_deref() {
        let d = cmp::min(record.delta, 255) as usize;
        ctx += u32::from(d_tab[d]) << param.d_loc;

        if record.prevq != q {
            record.delta += 1;
        }

        record.prevq = q;
    }

    if param.flags.contains(Flags::DO_SEL) {
        ctx += u32::from(record.sel) << param.s_loc;
    }

    (ctx & 0xffff) as u16
}

fn decode_length(
    src: &mut &[u8],
    range_coder: &mut RangeCoder,
    models: &mut Models,
) -> io::Result<u32> {
    let b0 = models.len[0].decode(src, range_coder).map(u32::from)?;
    let b1 = models.len[1].decode(src, range_coder).map(u32::from)?;
    let b2 = models.len[2].decode(src, range_coder).map(u32::from)?;
    let b3 = models.len[3].decode(src, range_coder).map(u32::from)?;

    Ok((b3 << 24) | (b2 << 16) | (b1 << 8) | b0)
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
}
