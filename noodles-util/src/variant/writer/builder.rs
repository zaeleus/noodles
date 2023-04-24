use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;

use crate::variant::{Compression, Format};

use super::Writer;

/// A variant writer builder.
#[derive(Default)]
pub struct Builder {
    // None means infer on build; Some(None) means no compression explicitly set.
    compression: Option<Option<Compression>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Compression};
    /// let builder = variant::writer::Builder::default().set_compression(Some(Compression::Bgzf));
    /// ```
    pub fn set_compression(mut self, compression: Option<Compression>) -> Self {
        self.compression = Some(compression);
        self
    }

    /// Sets the format of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, Format};
    /// let builder = variant::writer::Builder::default().set_format(Format::Vcf);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Builds a variant writer from a path.
    ///
    /// If the format or compression is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::variant::{self, Compression, Format};
    ///
    /// let writer = variant::writer::Builder::default()
    ///     .set_format(Format::Vcf)
    ///     .set_compression(Some(Compression::Bgzf))
    ///     .build_from_path("out.vcf.gz")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let compression = self
            .compression
            .get_or_insert_with(|| detect_compression_from_path_extension(src));

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src, *compression);
        }

        let file = File::create(src).map(BufWriter::new)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds a variant writer from a writer.
    ///
    /// If the format is not set, VCF format is used. If the compression is not set, no compression
    /// is used.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::variant::{self, Compression, Format};
    ///
    /// let writer = variant::writer::Builder::default()
    ///     .set_format(Format::Vcf)
    ///     .set_compression(Some(Compression::Bgzf))
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Vcf);
        let compression = self.compression.unwrap_or(match format {
            Format::Vcf => None,
            Format::Bcf => Some(Compression::Bgzf),
        });

        let inner: Box<dyn vcf::VariantWriter> = match (format, compression) {
            (Format::Vcf, None) => Box::new(vcf::Writer::new(writer)),
            (Format::Vcf, Some(Compression::Bgzf)) => {
                Box::new(vcf::Writer::new(bgzf::Writer::new(writer)))
            }
            (Format::Bcf, None) => Box::new(bcf::Writer::from(writer)),
            (Format::Bcf, Some(Compression::Bgzf)) => Box::new(bcf::Writer::new(writer)),
        };

        Writer { inner }
    }
}

fn detect_format_from_path_extension<P>(path: P, compression: Option<Compression>) -> Option<Format>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();
    let ext = path.extension().and_then(|ext| ext.to_str());

    match (compression, ext) {
        (None, Some("vcf")) => Some(Format::Vcf),
        (Some(Compression::Bgzf), Some("gz" | "bgz")) => {
            let path: &Path = path.file_stem()?.as_ref();
            let ext = path.extension().and_then(|ext| ext.to_str());

            match ext {
                Some("vcf") => Some(Format::Vcf),
                _ => None,
            }
        }
        (None | Some(Compression::Bgzf), Some("bcf")) => Some(Format::Bcf),
        _ => None,
    }
}

fn detect_compression_from_path_extension<P>(path: P) -> Option<Compression>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bcf" | "gz" | "bgz") => Some(Compression::Bgzf),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_from_path_extension() {
        assert_eq!(detect_compression_from_path_extension("out.vcf"), None);
        assert_eq!(
            detect_compression_from_path_extension("out.vcf.gz"),
            Some(Compression::Bgzf),
        );
        assert_eq!(
            detect_compression_from_path_extension("out.bcf"),
            Some(Compression::Bgzf)
        );
    }

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.vcf", None),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.vcf.gz", Some(Compression::Bgzf)),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bcf", None),
            Some(Format::Bcf)
        );

        assert_eq!(
            detect_format_from_path_extension("out.bcf.gz", Some(Compression::Bgzf)),
            None
        );
        assert_eq!(detect_format_from_path_extension("out.vcf.gz", None), None);

        assert!(detect_format_from_path_extension("out.fa", None).is_none());
    }
}
