mod flags;
pub mod parameter;

use std::{io, num::NonZero};

pub use self::{flags::Flags, parameter::Parameter};
use self::{flags::read_flags, parameter::fqz_decode_single_param};
use crate::io::reader::num::read_u8;

const VERSION: u8 = 5;

pub struct Parameters {
    pub gflags: Flags,
    pub max_sel: u8,
    pub s_tab: Vec<u8>,
    pub params: Vec<Parameter>,
    pub max_symbol_count: NonZero<usize>,
}

pub fn fqz_decode_params(src: &mut &[u8]) -> io::Result<Parameters> {
    let vers = read_u8(src)?;

    if vers != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid vers: expected {VERSION}, got {vers}"),
        ));
    }

    let gflags = read_flags(src)?;

    let (n_param, mut max_sel) = if gflags.contains(Flags::MULTI_PARAM) {
        let n = read_u8(src)?;
        (usize::from(n), n)
    } else {
        (1, 0)
    };

    let s_tab = if gflags.contains(Flags::HAVE_S_TAB) {
        max_sel = read_u8(src)?;
        read_array(src, 256)?
    } else {
        Vec::new()
    };

    let mut params = Vec::with_capacity(n_param);

    for _ in 0..n_param {
        let param = fqz_decode_single_param(src)?;
        params.push(param);
    }

    let max_symbol_count = params.iter().map(|param| param.symbol_count).max().unwrap();

    Ok(Parameters {
        gflags,
        max_sel,
        s_tab,
        params,
        max_symbol_count,
    })
}

pub fn read_array(src: &mut &[u8], n: usize) -> io::Result<Vec<u8>> {
    let (mut j, mut z) = (0, 0);
    let mut last = 0;

    let mut runs = vec![0; n];

    while z < n {
        let run = read_u8(src)?;

        runs[j] = run;
        j += 1;
        z += usize::from(run);

        if run == last {
            let copy = read_u8(src)?;

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
