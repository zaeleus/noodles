mod flags;

use std::io::{self, Read};

use byteorder::ReadBytesExt;

use self::flags::Flags;
use super::{
    parameter::{fqz_decode_single_param, Parameter},
    read_array,
};

const VERSION: u8 = 5;

pub struct Parameters {
    gflags: Flags,
    pub max_sel: u8,
    s_tab: Vec<u8>,
    params: Vec<Parameter>,
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
            format!("invalid vers: expected {}, got {}", VERSION, vers),
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
