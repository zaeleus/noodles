mod flags;

pub use self::flags::Flags;

use std::io::{self, Read};

use byteorder::ReadBytesExt;

use super::parameter::{Parameter, fqz_decode_single_param};

const VERSION: u8 = 5;

pub struct Parameters {
    pub gflags: Flags,
    pub max_sel: u8,
    pub s_tab: Vec<u8>,
    pub params: Vec<Parameter>,
    pub max_sym: u8,
}

pub fn fqz_decode_params<R>(reader: &mut R) -> io::Result<Parameters>
where
    R: Read,
{
    let vers = reader.read_u8()?;

    if vers != VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid vers: expected {VERSION}, got {vers}"),
        ));
    }

    let gflags = reader.read_u8().map(Flags::from)?;

    let (n_param, mut max_sel) = if gflags.contains(Flags::MULTI_PARAM) {
        let n = reader.read_u8()?;
        (usize::from(n), n)
    } else {
        (1, 0)
    };

    let s_tab = if gflags.contains(Flags::HAVE_S_TAB) {
        max_sel = reader.read_u8()?;
        read_array(reader, 256)?
    } else {
        Vec::new()
    };

    let mut params = Vec::with_capacity(n_param);
    let mut max_sym = 0;

    for _ in 0..n_param {
        let param = fqz_decode_single_param(reader)?;

        if param.max_sym > max_sym {
            max_sym = param.max_sym;
        }

        params.push(param);
    }

    Ok(Parameters {
        gflags,
        max_sel,
        s_tab,
        params,
        max_sym,
    })
}

pub fn read_array<R>(reader: &mut R, n: usize) -> io::Result<Vec<u8>>
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
