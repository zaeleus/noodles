use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use noodles_bcf as bcf;
use noodles_bgzf as bgzf;
use noodles_vcf as vcf;

use super::Writer;
#[allow(deprecated)]
use crate::variant::Compression;
use crate::variant::{CompressionMethod, Format};

/// A variant writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
}

impl Builder {
    /// Sets the compression method of the output.
    #[allow(deprecated)]
    #[deprecated(since = "0.20.0", note = "Use `Self::set_compression_method` instead.")]
    pub fn set_compression(self, compression: Option<Compression>) -> Self {
        self.set_compression_method(compression)
    }

    /// Sets the compression method of the output.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::variant::{self, CompressionMethod};
    /// let builder = variant::writer::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
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
    /// If the format or compression method is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::variant::{self, Format};
    /// let writer = variant::writer::Builder::default().build_from_path("out.vcf.gz")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        let compression = self
            .compression_method
            .get_or_insert_with(|| detect_compression_method_from_path_extension(src));

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src, *compression);
        }

        let file = File::create(src).map(BufWriter::new)?;
        Ok(self.build_from_writer(file))
    }

    /// Builds a variant writer from a writer.
    ///
    /// If the format is not set, VCF format is used. If the compression method is not set, a
    /// default one is used depending on the format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::variant::{self, Compression, Format};
    ///
    /// let writer = variant::writer::Builder::default()
    ///     .set_format(Format::Vcf)
    ///     .build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> Writer
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Vcf);
        let compression = self.compression_method.unwrap_or(match format {
            Format::Vcf => None,
            Format::Bcf => Some(CompressionMethod::Bgzf),
        });

        let inner: Box<dyn vcf::VariantWriter> = match (format, compression) {
            (Format::Vcf, None) => Box::new(vcf::Writer::new(writer)),
            (Format::Vcf, Some(CompressionMethod::Bgzf)) => {
                Box::new(vcf::Writer::new(bgzf::Writer::new(writer)))
            }
            (Format::Bcf, None) => Box::new(bcf::Writer::from(writer)),
            (Format::Bcf, Some(CompressionMethod::Bgzf)) => Box::new(bcf::Writer::new(writer)),
        };

        Writer { inner }
    }
}

fn detect_format_from_path_extension<P>(
    path: P,
    compression: Option<CompressionMethod>,
) -> Option<Format>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();
    let ext = path.extension().and_then(|ext| ext.to_str());

    match (compression, ext) {
        (None, Some("vcf")) => Some(Format::Vcf),
        (Some(CompressionMethod::Bgzf), Some("gz" | "bgz")) => {
            let path: &Path = path.file_stem()?.as_ref();
            let ext = path.extension().and_then(|ext| ext.to_str());

            match ext {
                Some("vcf") => Some(Format::Vcf),
                _ => None,
            }
        }
        (None | Some(CompressionMethod::Bgzf), Some("bcf")) => Some(Format::Bcf),
        _ => None,
    }
}

fn detect_compression_method_from_path_extension<P>(path: P) -> Option<CompressionMethod>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bcf" | "gz" | "bgz") => Some(CompressionMethod::Bgzf),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_method_from_path_extension() {
        assert_eq!(
            detect_compression_method_from_path_extension("out.vcf"),
            None
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.vcf.gz"),
            Some(CompressionMethod::Bgzf),
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.bcf"),
            Some(CompressionMethod::Bgzf)
        );
    }

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.vcf", None),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.vcf.gz", Some(CompressionMethod::Bgzf)),
            Some(Format::Vcf)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bcf", None),
            Some(Format::Bcf)
        );

        assert_eq!(
            detect_format_from_path_extension("out.bcf.gz", Some(CompressionMethod::Bgzf)),
            None
        );
        assert_eq!(detect_format_from_path_extension("out.vcf.gz", None), None);

        assert!(detect_format_from_path_extension("out.fa", None).is_none());
    }
}
