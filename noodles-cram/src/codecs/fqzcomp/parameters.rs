mod flags;
pub mod parameter;

use std::{io, num::NonZero};

pub use self::{flags::Flags, parameter::Parameter};
use self::{flags::read_flags, parameter::fqz_decode_single_param};
use crate::io::reader::num::read_u8;

pub struct Parameters {
    pub gflags: Flags,
    pub max_sel: u8,
    pub s_tab: Vec<u8>,
    pub params: Vec<Parameter>,
    pub max_symbol_count: NonZero<usize>,
}

pub fn fqz_decode_params(src: &mut &[u8]) -> io::Result<Parameters> {
    read_version(src)?;

    let gflags = read_flags(src)?;

    let (parameter_count, mut max_sel) = if gflags.contains(Flags::MULTI_PARAM) {
        let n = read_parameter_count(src)?;
        let selector_count = (usize::from(n) - 1) as u8;
        (n, selector_count)
    } else {
        (NonZero::<usize>::MIN, 0)
    };

    let s_tab = if gflags.contains(Flags::HAVE_S_TAB) {
        max_sel = read_u8(src)?;
        read_array(src, 256)?
    } else {
        Vec::new()
    };

    let params: Vec<_> = (0..parameter_count.get())
        .map(|_| fqz_decode_single_param(src))
        .collect::<io::Result<_>>()?;

    // SAFETY: `params` is nonempty.
    let max_symbol_count = params.iter().map(|param| param.symbol_count).max().unwrap();

    Ok(Parameters {
        gflags,
        max_sel,
        s_tab,
        params,
        max_symbol_count,
    })
}

fn read_version(src: &mut &[u8]) -> io::Result<()> {
    const VERSION: u8 = 5;

    match read_u8(src)? {
        VERSION => Ok(()),
        n => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid version: expected {VERSION}, got {n}"),
        )),
    }
}

fn read_parameter_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_u8(src).and_then(|n| {
        NonZero::try_from(usize::from(n)).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
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
