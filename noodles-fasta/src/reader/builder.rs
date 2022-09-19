use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufRead, BufReader},
    path::{Path, PathBuf},
};

use noodles_bgzf::{self as bgzf, gzi};

use crate::{fai, reader::Source, Reader};

/// A FASTA reader builder.
#[derive(Debug, Default)]
pub struct Builder {
    fasta_index: Option<fai::Index>,
    gz_index: Option<gzi::Index>,
}

impl Builder {
    /// Set the FASTA index.
    pub fn set_fasta_index(mut self, index: fai::Index) -> Self {
        self.fasta_index = Some(index);
        self
    }

    /// Sets the FASTA bgzip index.
    pub fn set_gz_index(mut self, index: gzi::Index) -> Self {
        self.gz_index = Some(index);
        self
    }

    /// Builds a FASTA reader from a path.
    ///
    /// If an associated FASTA index exists, it will be read.
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Reader<Source<File>>>
    where
        P: AsRef<Path>,
    {
        if self.fasta_index.is_none() {
            let fasta_index_src = build_fasta_index_src(&src);

            if fasta_index_src.try_exists()? {
                let index = fai::read(fasta_index_src)?;
                self = self.set_fasta_index(index);
            }
        }

        let file = File::open(&src)?;
        let source = match src.as_ref().extension().and_then(|ext| ext.to_str()) {
            Some("gz") => {
                if self.gz_index.is_none() {
                    let gz_index_src = build_gz_index_src(src);

                    if gz_index_src.try_exists()? {
                        let index = gzi::read(gz_index_src)?;
                        self = self.set_gz_index(index);
                    }
                }

                let reader = bgzf::Reader::new(file);

                match self.gz_index.take() {
                    Some(index) => Source::IndexedBgzip(reader, index),
                    _ => Source::from(reader),
                }
            }
            _ => Source::from(BufReader::new(file)),
        };

        self.build_from_reader(source)
    }

    /// Builds a FASTA reader from a reader.
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<R>>
    where
        R: BufRead,
    {
        use super::inner::{IndexedRawReader, Inner, RawReader};

        let inner = if let Some(index) = self.fasta_index {
            Inner::IndexedRaw(IndexedRawReader::new(reader, index))
        } else {
            Inner::Raw(RawReader::new(reader))
        };

        Ok(Reader { inner })
    }
}

fn build_fasta_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "fai";
    push_ext(src.as_ref().into(), EXT)
}

fn build_gz_index_src<P>(src: P) -> PathBuf
where
    P: AsRef<Path>,
{
    const EXT: &str = "gzi";
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
    fn test_build_fasta_index_src() {
        assert_eq!(build_fasta_index_src("ref.fa"), PathBuf::from("ref.fa.fai"));
    }

    #[test]
    fn test_build_gz_index_src() {
        assert_eq!(
            build_gz_index_src("ref.fa.gz"),
            PathBuf::from("ref.fa.gz.gzi")
        );
    }
}
