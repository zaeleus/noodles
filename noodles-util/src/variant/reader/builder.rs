use std::{
    fs::File,
    io::{self, BufRead, BufReader, Read},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf::{self as vcf, VariantReader};

use crate::variant::{Compression, Format};

use super::Reader;

/// A variant reader builder.
#[derive(Default)]
pub struct Builder {
    // None means infer on build; Some(None) means no compression explicitly set.
    compression: Option<Option<Compression>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression of the input.
    ///
    /// By default, the compression is autodetected on build. This can be used to override it.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Compression};
    /// let builder = variant::reader::Builder::default().set_compression(Some(Compression::Bgzf));
    /// ```
    pub fn set_compression(mut self, compression: Option<Compression>) -> Self {
        self.compression = Some(compression);
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
    /// By default, the format and compression will be autodetected. This can be overridden by
    /// using [`Self::set_format`] and [`Self::set_compression`].
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
    /// By default, the format and compression will be autodetected. This can be overridden by
    /// using [`Self::set_format`] and [`Self::set_compression`].
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

        let compression = self
            .compression
            .map(Ok)
            .unwrap_or_else(|| detect_compression(&mut reader))?;

        let format = self
            .format
            .map(Ok)
            .unwrap_or_else(|| detect_format(&mut reader, compression))?;

        let inner: Box<dyn VariantReader<_>> = match (format, compression) {
            (Format::Vcf, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(vcf::Reader::new(inner))
            }
            (Format::Vcf, Some(Compression::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(vcf::Reader::new(inner))
            }
            (Format::Bcf, None) => {
                let inner: Box<dyn BufRead> = Box::new(reader);
                Box::new(bcf::Reader::from(inner))
            }
            (Format::Bcf, Some(Compression::Bgzf)) => {
                let inner: Box<dyn BufRead> = Box::new(bgzf::Reader::new(reader));
                Box::new(bcf::Reader::from(inner))
            }
        };

        Ok(Reader { inner })
    }
}

fn detect_compression<R>(reader: &mut R) -> io::Result<Option<Compression>>
where
    R: BufRead,
{
    const GZIP_MAGIC_NUMBER: [u8; 2] = [0x1f, 0x8b];

    let src = reader.fill_buf()?;
    if src.get(..2) == Some(&GZIP_MAGIC_NUMBER) {
        Ok(Some(Compression::Bgzf))
    } else {
        Ok(None)
    }
}

fn detect_format<R>(reader: &mut R, compression: Option<Compression>) -> io::Result<Format>
where
    R: BufRead,
{
    const BCF_MAGIC_NUMBER: [u8; 3] = *b"BCF";

    let src = reader.fill_buf()?;

    if compression == Some(Compression::Bgzf) {
        let mut reader = bgzf::Reader::new(src);
        let mut buf = [0; 3];
        reader.read_exact(&mut buf).ok();

        if buf == BCF_MAGIC_NUMBER {
            return Ok(Format::Bcf);
        }
    } else if src.get(..3) == Some(&BCF_MAGIC_NUMBER) {
        return Ok(Format::Bcf);
    }

    Ok(Format::Vcf)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_format() -> io::Result<()> {
        use std::io::Write;

        fn t(mut src: &[u8], compression: Option<Compression>, expected: Format) {
            assert!(matches!(detect_format(&mut src, compression), Ok(value) if value == expected));
        }

        let header = vcf::Header::default();
        let raw_header = header.to_string();

        let src = raw_header.as_bytes();
        t(src, None, Format::Vcf);

        let mut writer = bgzf::Writer::new(Vec::new());
        writer.write_all(raw_header.as_bytes())?;
        let src = writer.finish()?;
        t(&src, Some(Compression::Bgzf), Format::Vcf);

        let mut writer = bcf::Writer::from(Vec::new());
        writer.write_header(&header)?;
        let src = writer.into_inner();
        t(&src, None, Format::Bcf);

        let mut writer = bcf::Writer::new(Vec::new());
        writer.write_header(&header)?;
        let src = writer.into_inner().finish()?;
        t(&src, Some(Compression::Bgzf), Format::Bcf);

        Ok(())
    }
}
