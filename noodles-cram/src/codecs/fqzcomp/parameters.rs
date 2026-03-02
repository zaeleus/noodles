mod flags;
pub mod parameter;

use std::{io, num::NonZero};

pub use self::{flags::Flags, parameter::Parameter};
use self::{flags::read_flags, parameter::fqz_decode_single_param};
use crate::io::reader::num::read_u8;

// § 6.2 "FQZComp Data Stream" (2023-03-15).
const SELECTOR_TABLE_SIZE: usize = 256;

pub struct Parameters {
    pub gflags: Flags,
    selector_table: Option<Vec<u8>>,
    pub params: Vec<Parameter>,
    pub max_symbol_count: NonZero<usize>,
    pub selector_count: Option<NonZero<usize>>,
}

impl Parameters {
    pub(super) fn selector_table(&self) -> Option<&[u8]> {
        self.selector_table.as_deref()
    }
}

pub fn fqz_decode_params(src: &mut &[u8]) -> io::Result<Parameters> {
    read_version(src)?;

    let flags = read_flags(src)?;

    let mut selector_count = None;

    let parameter_count = if flags.has_parameter_count() {
        let n = read_parameter_count(src)?;
        selector_count = Some(n.checked_add(1).expect("attempt to add with overflow"));
        n
    } else {
        NonZero::<usize>::MIN
    };

    let selector_table = if flags.has_selector_table() {
        selector_count = read_selector_count(src).map(Some)?;
        read_selector_table(src).map(Some)?
    } else {
        None
    };

    let params: Vec<_> = (0..parameter_count.get())
        .map(|_| fqz_decode_single_param(src))
        .collect::<io::Result<_>>()?;

    // SAFETY: `params` is nonempty.
    let max_symbol_count = params
        .iter()
        .map(|param| param.symbol_count())
        .max()
        .unwrap();

    Ok(Parameters {
        gflags: flags,
        selector_table,
        params,
        max_symbol_count,
        selector_count,
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

fn read_max_selector(src: &mut &[u8]) -> io::Result<NonZero<u8>> {
    read_u8(src).and_then(|n| {
        NonZero::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    })
}

fn read_selector_count(src: &mut &[u8]) -> io::Result<NonZero<usize>> {
    read_max_selector(src).map(|n| {
        NonZero::<usize>::from(n)
            .checked_add(1)
            .expect("attempt to add with overflow")
    })
}

fn read_selector_table(src: &mut &[u8]) -> io::Result<Vec<u8>> {
    read_array(src, SELECTOR_TABLE_SIZE)
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

    let mut i: u16 = 0;
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
            let value = u8::try_from(i).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("run-length array index overflow: i = {i}"),
                )
            })?;
            a[z] = value;
            z += 1;
        }

        i += 1;
    }

    Ok(a)
}
