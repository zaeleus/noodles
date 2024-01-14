use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Reader;
use crate::alignment::{CompressionMethod, Format};

/// An alignment reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// By default, the compression method is autodetected on build. This can be used to override
    /// it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, CompressionMethod};
    /// let builder = alignment::reader::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::reader::Builder::default().set_format(Format::Sam);
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
    /// let builder = alignment::reader::Builder::default()
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
    /// [`Self::set_format`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::reader::Builder::default().build_from_path("sample.bam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, path: P) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        self.build_from_reader(file)
    }

    /// Builds an alignment reader from a reader.
    ///
    /// By default, the format will be autodetected. This can be overridden by using
    /// [`Self::set_format`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment;
    /// let reader = alignment::reader::Builder::default().build_from_reader(io::empty())?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_reader<R>(self, reader: R) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        R: Read + 'static,
    {
        let mut reader = BufReader::new(reader);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => detect_compression_method(&mut reader)?,
        };

        let format = match self.format {
            Some(format) => format,
            None => detect_format(&mut reader, compression_method)?,
        };

        let inner: Box<dyn sam::alignment::io::Read<_>> = match (format, compression_method) {
            (Format::Sam, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(sam::io::Reader::from(inner))
            }
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(sam::io::Reader::from(inner))
            }
            (Format::Bam, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(bam::io::Reader::from(inner))
            }
            (Format::Bam, Some(CompressionMethod::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(bam::io::Reader::from(inner))
            }
            (Format::Cram, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);

                Box::new(
                    cram::io::reader::Builder::default()
                        .set_reference_sequence_repository(self.reference_sequence_repository)
                        .build_from_reader(inner),
                )
            }
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "CRAM cannot be bgzip-compressed",
                ))
            }
        };

        Ok(Reader { inner })
    }
}

pub(crate) fn detect_compression_method<R>(reader: &mut R) -> io::Result<Option<CompressionMethod>>
where
    R: BufRead,
{
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let src = reader.fill_buf()?;

    if let Some(buf) = src.get(..GZIP_MAGIC_NUMBER.len()) {
        if buf == GZIP_MAGIC_NUMBER {
            return Ok(Some(CompressionMethod::Bgzf));
        }
    }

    Ok(None)
}

pub(crate) fn detect_format<R>(
    reader: &mut R,
    compression_method: Option<CompressionMethod>,
) -> io::Result<Format>
where
    R: BufRead,
{
    use flate2::bufread::MultiGzDecoder;

    const CRAM_MAGIC_NUMBER: [u8; 4] = [b'C', b'R', b'A', b'M'];
    const BAM_MAGIC_NUMBER: [u8; 4] = [b'B', b'A', b'M', 0x01];

    let src = reader.fill_buf()?;

    if matches!(compression_method, Some(CompressionMethod::Bgzf)) {
        let mut decoder = MultiGzDecoder::new(src);
        let mut buf = [0; BAM_MAGIC_NUMBER.len()];
        decoder.read_exact(&mut buf)?;

        if buf == BAM_MAGIC_NUMBER {
            return Ok(Format::Bam);
        }
    } else if let Some(buf) = src.get(..BAM_MAGIC_NUMBER.len()) {
        if buf == BAM_MAGIC_NUMBER {
            return Ok(Format::Bam);
        } else if buf == CRAM_MAGIC_NUMBER {
            return Ok(Format::Cram);
        }
    }

    Ok(Format::Sam)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_format() -> io::Result<()> {
        use std::io::Write;

        fn t(mut src: &[u8], expected: Format, compression_method: Option<CompressionMethod>) {
            assert!(
                matches!(detect_format(&mut src, compression_method), Ok(value) if value == expected)
            );
        }

        t(b"@HD\tVN:1.6\n", Format::Sam, None);
        t(b"", Format::Sam, None);

        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(b"@HD\tVN:1.6\n")?;
        let src = writer.finish()?;
        t(&src, Format::Sam, Some(CompressionMethod::Bgzf));

        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(b"BAM\x01")?;
        let src = dbg!(writer.finish())?;
        t(&src, Format::Bam, Some(CompressionMethod::Bgzf));

        // An incomplete gzip stream. See #179.
        #[rustfmt::skip]
        let src = [
            0x1f, 0x8b, // ID1, ID2
            0x08, // CM = DEFLATE
            0x04, // FLG = FEXTRA
            0x00, 0x00, 0x00, 0x00, // MTIME = 0
            0x00, // XFL = 0
            0xff, // OS = 255 (unknown)
            0x06, 0x00, // XLEN = 6
            b'B', b'C', // SI1, SI2
            0x02, 0x00, // SLEN = 2
            0x00, 0x40, // BSIZE = 16384
            0x73, 0x72, 0xf4, 0x65, 0x04, 0x00, // CDATA = deflate(b"BAM\x01")
            // ...
        ];
        t(&src, Format::Bam, Some(CompressionMethod::Bgzf));

        t(b"CRAM", Format::Cram, None);

        Ok(())
    }
}
