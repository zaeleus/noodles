use std::io::{self, Write};

use super::{
    Models,
    parameters::{self, parameter},
};
use crate::{
    codecs::aac::RangeCoder,
    io::writer::num::{write_u8, write_u16_le, write_uint7},
};

pub fn encode(lens: &[usize], src: &[u8]) -> io::Result<Vec<u8>> {
    let mut dst = Vec::new();

    let len =
        u32::try_from(src.len()).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    write_uint7(&mut dst, len)?;

    let mut q_hist = [[0; 256]];

    for (i, f) in q_hist[0].iter_mut().enumerate() {
        *f = i as u8;
    }

    let parameters = build_parameters(lens, src);
    fqz_encode_params(&mut dst, &parameters)?;

    let mut range_coder = RangeCoder::default();
    let mut models = Models::new(parameters.max_sym, Some(1)); // FIXME: max_sel

    let mut p = 0;
    let mut rec_num = 0;

    let mut x = 0;
    let mut last = 0;
    let mut qlast: u32 = 0;

    for &q in src {
        if p == 0 {
            if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
                todo!("have_s_tab");
            }

            x = parameters.s_tab[0];

            let is_fixed_len = parameters.params[usize::from(x)]
                .flags
                .contains(parameter::Flags::DO_LEN);

            let len = lens[rec_num];

            if !is_fixed_len || rec_num == 0 {
                encode_length(&mut dst, &mut range_coder, &mut models, len)?;
            }

            let param = &parameters.params[usize::from(x)];

            if param.flags.contains(parameter::Flags::DO_DEDUP) {
                todo!("do_dedup");
            }

            p = len;
            last = u32::from(param.context);
            qlast = 0;

            rec_num += 1;
        }

        let qq = q_hist[usize::from(x)][usize::from(q)];
        models.qual[usize::from(last as u16)].encode(&mut dst, &mut range_coder, qq)?;

        let param = &parameters.params[usize::from(x)];

        qlast = (qlast << param.q_shift)
            .overflowing_add(u32::from(param.q_tab[usize::from(qq)]))
            .0;
        last = u32::from(param.context);
        last += (qlast & ((1 << param.q_bits) - 1)) << param.q_loc;

        if param.flags.contains(parameter::Flags::HAVE_PTAB) {
            last += u32::from(param.p_tab[p.min(1023)] << param.p_loc);
        }

        if param.flags.contains(parameter::Flags::HAVE_DTAB) {
            todo!("have_dtab");
        }

        if param.flags.contains(parameter::Flags::DO_SEL) {
            todo!("do_sel");
        }

        last &= 0xffff;
        p -= 1;
    }

    range_coder.range_encode_end(&mut dst)?;

    Ok(dst)
}

struct Parameters {
    pub gflags: parameters::Flags,
    pub max_sel: u8,
    pub s_tab: Vec<u8>,
    pub params: Vec<Parameter>,
    pub max_sym: u8,
}

struct Parameter {
    pub context: u16,
    pub flags: parameter::Flags,

    pub max_sym: u8,

    pub q_bits: u8,
    pub q_shift: u8,
    pub q_loc: u8,

    pub s_loc: u8,

    pub p_loc: u8,

    pub d_loc: u8,

    pub q_tab: Vec<u8>,
    pub p_tab: Vec<u8>,
}

fn build_parameters(lens: &[usize], src: &[u8]) -> Parameters {
    let mut max_symbol = u8::MIN;
    let mut symbol_counts = [0; 256];

    for &b in src {
        let sym = usize::from(b);
        symbol_counts[sym] += 1;
        max_symbol = max_symbol.max(b);
    }

    let q_shift = 5;
    let q_bits = if q_shift > 4 { 9 } else { 8 };
    let p_bits = 7;
    let p_shift = i32::from(lens[0] > 128);

    let q_tab: Vec<_> = (0..=u8::MAX).collect();

    let mut p_tab = vec![0; 1024];

    for (i, p) in p_tab.iter_mut().enumerate() {
        *p = ((1 << p_bits) - 1).min(i >> p_shift) as u8;
    }

    let mut flags = parameter::Flags::HAVE_PTAB;

    if lens.windows(2).all(|w| w[0] == w[1]) {
        flags |= parameter::Flags::DO_LEN;
    }

    let params = vec![Parameter {
        context: 0,
        flags,
        max_sym: max_symbol,
        q_bits,
        q_shift,
        q_loc: 7,
        s_loc: 15,
        p_loc: 0,
        d_loc: 15,
        q_tab,
        p_tab,
    }];

    let last_i = (params.len() - 1) as u8;
    let mut s_tab = vec![last_i; 256];

    for (i, s) in s_tab.iter_mut().enumerate().take(params.len()) {
        *s = i as u8;
    }

    Parameters {
        gflags: parameters::Flags::empty(),
        max_sel: 0,
        s_tab,
        params,
        max_sym: max_symbol,
    }
}

fn fqz_encode_params<W>(writer: &mut W, parameters: &Parameters) -> io::Result<()>
where
    W: Write,
{
    const VERSION: u8 = 5;

    write_u8(writer, VERSION)?;

    let gflags = u8::from(parameters.gflags);
    write_u8(writer, gflags)?;

    if parameters.gflags.contains(parameters::Flags::MULTI_PARAM) {
        let n_param = u8::try_from(parameters.params.len())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        write_u8(writer, n_param)?;
    }

    if parameters.gflags.contains(parameters::Flags::HAVE_S_TAB) {
        write_u8(writer, parameters.max_sel)?;
        write_array(writer, &parameters.s_tab)?;
    }

    for param in &parameters.params {
        fqz_encode_single_param(writer, param)?;
    }

    Ok(())
}

fn fqz_encode_single_param<W>(writer: &mut W, parameter: &Parameter) -> io::Result<()>
where
    W: Write,
{
    let context = parameter.context;
    write_u16_le(writer, context)?;

    let pflags = u8::from(parameter.flags);
    write_u8(writer, pflags)?;

    let max_sym = parameter.max_sym;
    write_u8(writer, max_sym)?;

    write_u8(writer, (parameter.q_bits << 4) | parameter.q_shift)?;
    write_u8(writer, (parameter.q_loc << 4) | parameter.s_loc)?;
    write_u8(writer, (parameter.p_loc << 4) | parameter.d_loc)?;

    if parameter.flags.contains(parameter::Flags::HAVE_QMAP) {
        todo!("have_qmap");
    }

    if parameter.flags.contains(parameter::Flags::HAVE_QTAB) {
        todo!("have_qtab");
    }

    if parameter.flags.contains(parameter::Flags::HAVE_PTAB) {
        write_array(writer, &parameter.p_tab)?;
    }

    if parameter.flags.contains(parameter::Flags::HAVE_DTAB) {
        todo!("have_dtab");
    }

    Ok(())
}

fn write_array<W>(writer: &mut W, data: &[u8]) -> io::Result<()>
where
    W: Write,
{
    let mut rle1 = Vec::new();

    let mut i = 0;
    let mut j = 0;

    while j < data.len() {
        let start = j;

        while j < data.len() && usize::from(data[j]) == i {
            j += 1;
        }

        let mut len = j - start;

        loop {
            let rle = len.min(255);
            rle1.push(rle as u8);

            len -= rle;

            if rle != 255 {
                break;
            }
        }

        i += 1;
    }

    let mut rle2 = Vec::new();
    j = 0;
    let mut last = -1;

    while j < rle1.len() {
        let curr = rle1[j];
        j += 1;

        rle2.push(curr);

        if i32::from(curr) == last {
            let start = j;
            let mut len = 0;

            while j < rle1.len() && i32::from(rle1[j]) == last && len < 255 {
                j += 1;
                len = j - start;
            }

            rle2.push(len as u8);
        } else {
            last = i32::from(curr);
        }
    }

    writer.write_all(&rle2)?;

    Ok(())
}

fn encode_length<W>(
    writer: &mut W,
    range_coder: &mut RangeCoder,
    models: &mut Models,
    len: usize,
) -> io::Result<()>
where
    W: Write,
{
    let n = u32::try_from(len).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    models.len[0].encode(writer, range_coder, (n & 0xff) as u8)?;
    models.len[1].encode(writer, range_coder, ((n >> 8) & 0xff) as u8)?;
    models.len[2].encode(writer, range_coder, ((n >> 16) & 0xff) as u8)?;
    models.len[3].encode(writer, range_coder, ((n >> 24) & 0xff) as u8)?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0],
        ];

        let lens: Vec<_> = data.iter().map(|scores| scores.len()).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&lens, &src)?;

        let expected = [
            0x19, 0x05, 0x00, 0x00, 0x00, 0x20, 0x03, 0x95, 0x7f, 0x0f, 0x01, 0x01, 0x7d, 0xff,
            0xff, 0x01, 0x84, 0x00, 0x09, 0xff, 0xff, 0xf6, 0x01, 0x65, 0x00, 0x86, 0x2e, 0x98,
            0xea, 0xca, 0x71, 0x6f, 0x08, 0x51, 0x6f, 0x00,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_encode_with_do_len() -> io::Result<()> {
        let data = [
            vec![0, 0, 0, 1, 1, 2, 1, 1, 0, 0],
            vec![0, 1, 2, 3, 3, 3, 3, 3, 3, 3],
            vec![2, 1, 1, 0, 0, 0, 0, 0, 1, 1],
        ];

        let lens: Vec<_> = data.iter().map(|scores| scores.len()).collect();
        let src: Vec<_> = data.into_iter().flatten().collect();

        let actual = encode(&lens, &src)?;

        let expected = [
            0x1e, 0x05, 0x00, 0x00, 0x00, 0x24, 0x03, 0x95, 0x7f, 0x0f, 0x01, 0x01, 0x7d, 0xff,
            0xff, 0x01, 0x84, 0x00, 0x09, 0xff, 0xff, 0xf6, 0x01, 0x65, 0x0c, 0x10, 0x86, 0x6d,
            0x57, 0x10, 0x38, 0x60, 0xac,
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
