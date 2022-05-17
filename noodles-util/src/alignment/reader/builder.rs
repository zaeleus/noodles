use std::{
    ffi::{OsStr, OsString},
    fs::File,
    io::{self, BufReader, Read, Seek},
    path::{Path, PathBuf},
};

use noodles_bam::{self as bam, bai};
use noodles_bgzf as bgzf;
use noodles_cram::{self as cram, crai};
use noodles_csi as csi;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Reader;
use crate::alignment::Format;

/// An alignment reader builder.
pub struct Builder {
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    index_src: Option<PathBuf>,
}

impl Builder {
    pub(super) fn new() -> Self {
        Self {
            reference_sequence_repository: fasta::Repository::default(),
            format: None,
            index_src: None,
        }
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on [`build`]. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::Reader::builder().set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::Reader::builder()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Builds an alignment reader from a path.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`set_format`]. An associated index will also attempt to be loaded.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::Reader::builder().build_from_path("sample.bam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, path: P) -> io::Result<Reader<File>>
    where
        P: AsRef<Path>,
    {
        self.index_src = find_index_src(&path);
        let file = File::open(path)?;
        self.build_from_reader(file)
    }

    /// Builds an alignment reader from a reader.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`set_format`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::Reader::builder().build_from_reader(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, mut reader: R) -> io::Result<Reader<R>>
    where
        R: Read + Seek,
    {
        use super::{Index, Inner};

        let format = self
            .format
            .map(Ok)
            .unwrap_or_else(|| detect_format(&mut reader))?;

        let inner = match format {
            Format::Sam => Inner::Sam(sam::Reader::new(BufReader::new(reader))),
            Format::Bam => Inner::Bam(bam::Reader::new(reader)),
            Format::Cram => Inner::Cram(cram::Reader::new(reader)),
        };

        let mut index = None;

        if let Some(index_src) = self.index_src {
            index = match index_src.extension().and_then(|ext| ext.to_str()) {
                Some("bai") => bai::read(index_src).map(Index::Bai).map(Some)?,
                Some("crai") => crai::read(index_src).map(Index::Crai).map(Some)?,
                Some("csi") => csi::read(index_src).map(Index::Csi).map(Some)?,
                _ => None,
            }
        }

        Ok(Reader {
            inner,
            reference_sequence_repository: self.reference_sequence_repository,
            index,
        })
    }
}

fn detect_format<R>(reader: &mut R) -> io::Result<Format>
where
    R: Read + Seek,
{
    const CRAM_MAGIC_NUMBER: [u8; 4] = [b'C', b'R', b'A', b'M'];
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];
    const BAM_MAGIC_NUMBER: [u8; 4] = [b'B', b'A', b'M', 0x01];

    let mut buf = [0; 4];
    reader.read_exact(&mut buf).ok();
    reader.rewind()?;

    if buf == CRAM_MAGIC_NUMBER {
        return Ok(Format::Cram);
    }

    if buf[..2] == GZIP_MAGIC_NUMBER {
        let mut reader = bgzf::Reader::new(reader);
        reader.read_exact(&mut buf).ok();
        reader.get_mut().rewind()?;

        if buf == BAM_MAGIC_NUMBER {
            return Ok(Format::Bam);
        }
    }

    Ok(Format::Sam)
}

fn find_index_src<P>(src: P) -> Option<PathBuf>
where
    P: AsRef<Path>,
{
    const EXTENSIONS: [&str; 3] = ["bai", "crai", "csi"];

    let src = src.as_ref();

    for ext in EXTENSIONS {
        let index_src = push_ext(src.into(), ext);

        if index_src.exists() {
            return Some(index_src);
        }
    }

    None
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
