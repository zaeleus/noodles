#![allow(dead_code)]

mod parameter;
mod parameters;

use std::{
    cmp,
    io::{self, Read},
};

use byteorder::ReadBytesExt;

use self::{
    parameter::Parameter,
    parameters::{fqz_decode_params, Parameters},
};
use super::aac::{Model, RangeCoder};
use crate::reader::num::read_uint7;

pub fn fqz_decode<R>(reader: &mut R) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let buf_len = read_uint7(reader).and_then(|n| {
        usize::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })?;

    let mut params = fqz_decode_params(reader)?;

    let (mut range_coder, mut models) = fqz_create_models(&params);
    range_coder.range_decode_create(reader)?;

    let mut i = 0;
    let mut pos = 0;

    let mut record = Record::default();
    let mut dst = vec![0; buf_len];

    let mut x = 0;
    let mut ctx = 0;

    while i < buf_len {
        if pos == 0 {
            x = fqz_new_record(
                reader,
                &mut params,
                &mut range_coder,
                &mut models,
                &mut record,
            )?;

            if record.is_dup {
                for j in 0..record.rec_len {
                    dst[i + j] = dst[i + j - record.rec_len];
                }

                i += record.rec_len;
                pos = 0;

                continue;
            }

            ctx = params.params[x].context;
        }

        let param = &mut params.params[x];
        let q = models.qual[usize::from(ctx)].decode(reader, &mut range_coder)?;

        dst[i] = if param.flags.contains(parameter::Flags::HAVE_QMAP) {
            param.q_map[usize::from(q)]
        } else {
            q
        };

        ctx = fqz_update_context(param, q, &mut record);

        i += 1;
        pos -= 1;
    }

    if params.gflags.contains(parameters::Flags::DO_REV) {
        todo!("fqz_decode: reverse_qualities");
    }

    Ok(dst)
}

struct Models {
    len: Vec<Model>,
    qual: Vec<Model>,
    dup: Model,
    rev: Model,
    sel: Model,
}

fn fqz_create_models(parameters: &Parameters) -> (RangeCoder, Models) {
    let models = Models {
        len: vec![Model::new(u8::MAX); 4],
        qual: vec![Model::new(parameters.max_sym + 1); 1 << 16],
        dup: Model::new(2),
        rev: Model::new(2),
        sel: Model::new(parameters.max_sel + 1),
    };

    (RangeCoder::default(), models)
}

#[derive(Debug, Default)]
struct Record {
    rec: usize,
    sel: u8,
    rec_len: usize,
    is_dup: bool,
    qctx: u32,
    delta: u32,
    prevq: u8,
}

fn fqz_new_record<R>(
    reader: &mut R,
    parameters: &mut Parameters,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    record: &mut Record,
) -> io::Result<usize>
where
    R: Read,
{
    let mut sel = 0;
    let mut x = 0;

    if parameters.max_sel > 0 {
        sel = models.sel.decode(reader, range_coder)?;

        if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
            x = usize::from(parameters.s_tab[usize::from(sel)]);
        }
    }

    record.sel = sel;

    let param = &mut parameters.params[x];

    if param.flags.contains(parameter::Flags::DO_LEN) || param.first_len > 0 {
        let len = decode_length(reader, range_coder, models)?;
        param.last_len = len as usize;

        if !param.flags.contains(parameter::Flags::DO_LEN) {
            param.first_len = 0;
        }
    }

    record.rec_len = param.last_len;

    if parameters.gflags.contains(parameters::Flags::DO_REV) {
        todo!("fqz_new_record: do_rev");
    }

    record.rec += 1;

    if param.flags.contains(parameter::Flags::DO_DEDUP) {
        record.is_dup = models.dup.decode(reader, range_coder)? == 1;
    }

    record.qctx = 0;
    record.delta = 0;
    record.prevq = 0;

    Ok(x)
}

fn fqz_update_context(param: &mut Parameter, q: u8, record: &mut Record) -> u16 {
    use parameter::Flags;

    let mut ctx = u32::from(param.context);

    record.qctx =
        (record.qctx << u32::from(param.q_shift)) + u32::from(param.q_tab[usize::from(q)]);

    ctx += record.qctx & ((1 << param.q_bits) - 1) << param.q_loc;

    if param.flags.contains(Flags::HAVE_PTAB) {
        let p = cmp::min(record.rec_len, 1023) as usize;
        ctx += u32::from(param.p_tab[p]) << param.p_loc;
    }

    if param.flags.contains(Flags::HAVE_DTAB) {
        let d = cmp::min(record.delta, 255) as usize;
        ctx += u32::from(param.d_tab[d]) << param.d_loc;

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

fn decode_length<R>(
    reader: &mut R,
    range_coder: &mut RangeCoder,
    models: &mut Models,
) -> io::Result<u32>
where
    R: Read,
{
    let b0 = models.len[0].decode(reader, range_coder).map(u32::from)?;
    let b1 = models.len[1].decode(reader, range_coder).map(u32::from)?;
    let b2 = models.len[2].decode(reader, range_coder).map(u32::from)?;
    let b3 = models.len[3].decode(reader, range_coder).map(u32::from)?;

    Ok(b3 << 24 | b2 << 16 | b1 << 8 | b0)
}

fn read_array<R>(reader: &mut R, n: usize) -> io::Result<Vec<u8>>
where
    R: Read,
{
    let (mut j, mut z) = (0, 0);
    let mut last = 0;

    let mut runs = vec![0; n];

    while z < n {
        let run = reader.read_u8()?;

        runs[j] = run;
        j += 1;
        z += usize::from(run);

        if run == last {
            let copy = reader.read_u8()?;

            for _ in 0..copy {
                runs[j] = run;
                j += 1;
            }

            z += usize::from(run) * usize::from(copy);
        }

        last = run;
    }

    let mut a = vec![0; n];

    let mut i = 0;
    j = 0;
    z = 0;

    while z < n {
        let mut run_len = 0;

        loop {
            let part = runs[j];
            j += 1;
            run_len += usize::from(part);

            if part != 255 {
                break;
            }
        }

        for _ in 0..run_len {
            a[z] = i;
            z += 1;
        }

        i += 1;
    }

    Ok(a)
}
