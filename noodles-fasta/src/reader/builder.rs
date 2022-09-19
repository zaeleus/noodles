use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use noodles_bgzf as bgzf;

use crate::{fai, reader::Source, Reader};

/// A FASTA reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    index: Option<fai::Index>,
}

impl Builder {
    /// Set the FASTA index.
    pub fn set_index(mut self, index: fai::Index) -> Self {
        self.index = Some(index);
        self
    }

    /// Builds a FASTA reader from a path.
    ///
    /// If an associated FASTA index exists, it will be read.
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Reader<Source<File>>>
    where
        P: AsRef<Path>,
    {
        let index_src = build_index_src(&src);

        if index_src.try_exists()? {
            let index = fai::read(index_src)?;
            self = self.set_index(index);
        }

        let source = File::open(&src).map(|file| {
            match src.as_ref().extension().and_then(|ext| ext.to_str()) {
                Some("gz") => Source::from(bgzf::Reader::new(file)),
                _ => Source::from(BufReader::new(file)),
            }
        })?;

        self.build_from_reader(source)
    }

    /// Builds a FASTA reader from a reader.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<R>>
    where
        R: BufRead,
    {
        use super::inner::{IndexedRawReader, Inner, RawReader};

        let inner = if let Some(index) = self.index {
            Inner::IndexedRaw(IndexedRawReader::new(reader, index))
        } else {
            Inner::Raw(RawReader::new(reader))
        };

        Ok(Reader { inner })
    }
}

fn build_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "fai";
    push_ext(src.as_ref().into(), EXT)
}

fn push_ext<S>(path: PathBuf, ext: S) -> PathBuf
where
    S: AsRef<OsStr>,
{
    let mut s = OsString::from(path);
    s.push(".");
    s.push(ext);
    PathBuf::from(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_index_src() {
        assert_eq!(build_index_src("ref.fa"), PathBuf::from("ref.fa.fai"));
    }
}
