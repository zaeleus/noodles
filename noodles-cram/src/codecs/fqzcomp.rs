#![allow(dead_code)]

mod parameter;
mod parameters;

use std::io::{self, Read};

use byteorder::ReadBytesExt;

use self::parameters::Parameters;
use super::aac::{Model, RangeCoder};

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

struct Record {
    rec: usize,
    sel: u8,
    rec_len: usize,
    is_dup: bool,
    qctx: u32,
    delta: u32,
    prevq: u32,
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
            run_len += part;

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
