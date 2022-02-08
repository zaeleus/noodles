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
