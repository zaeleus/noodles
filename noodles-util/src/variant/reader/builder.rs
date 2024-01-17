use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, VariantReader};

use super::Reader;
use crate::variant::{CompressionMethod, Format};

/// A variant reader builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression method of the input.
    ///
    /// By default, the compression method is autodetected on build. This can be used to override
    /// it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, CompressionMethod};
    /// let builder = variant::reader::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression);
        self
    }

    /// Sets the format of the input.
    ///
    /// By default, the format is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Format};
    /// let builder = variant::reader::Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Builds a variant reader from a path.
    ///
    /// By default, the format and compression method will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::variant;
    /// let reader = variant::reader::Builder::default().build_from_path("sample.vcf")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(self, path: P) -> io::Result<Reader<Box<dyn BufRead>>>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        self.build_from_reader(file)
    }

    /// Builds a variant reader from a reader.
    ///
    /// By default, the format and compression methods will be autodetected. This can be overridden
    /// by using [`Self::set_format`] and [`Self::set_compression_method`].
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::variant;
    /// let reader = variant::reader::Builder::default().build_from_reader(io::empty())?;
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

        let inner: Box<dyn VariantReader<_>> = match (format, compression_method) {
            (Format::Vcf, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(vcf::Reader::new(inner))
            }
            (Format::Vcf, Some(CompressionMethod::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(vcf::Reader::new(inner))
            }
            (Format::Bcf, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(bcf::Reader::from(inner))
            }
            (Format::Bcf, Some(CompressionMethod::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(bcf::Reader::from(inner))
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

    const BCF_MAGIC_NUMBER: [u8; 3] = *b"BCF";

    let src = reader.fill_buf()?;

    if let Some(compression_method) = compression_method {
        if compression_method == CompressionMethod::Bgzf {
            let mut decoder = MultiGzDecoder::new(src);
            let mut buf = [0; BCF_MAGIC_NUMBER.len()];
            decoder.read_exact(&mut buf)?;

            if buf == BCF_MAGIC_NUMBER {
                return Ok(Format::Bcf);
            }
        }
    } else if let Some(buf) = src.get(..BCF_MAGIC_NUMBER.len()) {
        if buf == BCF_MAGIC_NUMBER {
            return Ok(Format::Bcf);
        }
    }

    Ok(Format::Vcf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_method() -> io::Result<()> {
        let mut src = &[0x1f, 0x8b][..];
        assert_eq!(
            detect_compression_method(&mut src)?,
            Some(CompressionMethod::Bgzf)
        );

        let mut src = &b"fileformat=VCFv4.4\n"[..];
        assert!(detect_compression_method(&mut src)?.is_none());

        let mut src = &[][..];
        assert!(detect_compression_method(&mut src)?.is_none());

        Ok(())
    }

    #[test]
    fn test_detect_format() -> io::Result<()> {
        use std::io::Write;

        fn t(mut src: &[u8], compression_method: Option<CompressionMethod>, expected: Format) {
            assert!(
                matches!(detect_format(&mut src, compression_method), Ok(value) if value == expected)
            );
        }

        let header = vcf::Header::default();
        let raw_header = header.to_string();

        let src = raw_header.as_bytes();
        t(src, None, Format::Vcf);

        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(raw_header.as_bytes())?;
        let src = writer.finish()?;
        t(&src, Some(CompressionMethod::Bgzf), Format::Vcf);

        let mut writer = bcf::Writer::from(Vec::new());
        writer.write_header(&header)?;
        let src = writer.into_inner();
        t(&src, None, Format::Bcf);

        let mut writer = bcf::Writer::new(Vec::new());
        writer.write_header(&header)?;
        let src = writer.into_inner().finish()?;
        t(&src, Some(CompressionMethod::Bgzf), Format::Bcf);

        // An incomplete gzip stream.
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
            0x73, 0x72, 0x76, 0x03, 0x00, // CDATA = deflate(b"BCF")
            // ...
        ];
        t(&src, Some(CompressionMethod::Bgzf), Format::Bcf);

        Ok(())
    }
}
