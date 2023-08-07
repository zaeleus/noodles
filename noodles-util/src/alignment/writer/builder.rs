//! Alignment writer builder.

use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use cram::data_container::BlockContentEncoderMap;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_cram as cram;
use noodles_fasta as fasta;
use noodles_sam as sam;

use super::Writer;
use crate::alignment::{CompressionMethod, Format};

/// An alignment writer builder.
#[derive(Default)]
pub struct Builder {
    compression_method: Option<Option<CompressionMethod>>,
    format: Option<Format>,
    reference_sequence_repository: fasta::Repository,
    block_content_encoder_map: BlockContentEncoderMap,
}

impl Builder {
    /// Sets the compression method.
    ///
    /// If not set, a default compression method is selected depending on the format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, CompressionMethod};
    /// let builder = alignment::writer::Builder::default()
    ///     .set_compression_method(Some(CompressionMethod::Bgzf));
    /// ```
    pub fn set_compression_method(mut self, compression_method: Option<CompressionMethod>) -> Self {
        self.compression_method = Some(compression_method);
        self
    }

    /// Sets the format of the output.
    ///
    /// If not set, a default format is used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_util::alignment::{self, Format};
    /// let builder = alignment::writer::Builder::default().set_format(Format::Sam);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = Some(format);
        self
    }

    /// Sets the reference sequence repository.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fasta as fasta;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let repository = fasta::Repository::default();
    ///
    /// let builder = alignment::writer::Builder::default()
    ///     .set_reference_sequence_repository(repository);
    /// ```
    pub fn set_reference_sequence_repository(
        mut self,
        reference_sequence_repository: fasta::Repository,
    ) -> Self {
        self.reference_sequence_repository = reference_sequence_repository;
        self
    }

    /// Sets the block content-encoder map.
    ///
    /// This is only used when the output format is CRAM.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_cram::data_container::BlockContentEncoderMap;
    /// use noodles_util::alignment::{self, Format};
    ///
    /// let builder = alignment::writer::Builder::default()
    ///     .set_block_content_encoder_map(BlockContentEncoderMap::default());
    /// ```
    pub fn set_block_content_encoder_map(
        mut self,
        block_content_encoder_map: BlockContentEncoderMap,
    ) -> Self {
        self.block_content_encoder_map = block_content_encoder_map;
        self
    }

    /// Builds an alignment writer from a path.
    ///
    /// If the format or compression method is not set, it is detected from the path extension.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    /// let writer = alignment::writer::Builder::default().build_from_path("out.sam")?;
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn build_from_path<P>(mut self, src: P) -> io::Result<Writer>
    where
        P: AsRef<Path>,
    {
        let src = src.as_ref();

        if self.compression_method.is_none() {
            self.compression_method = Some(detect_compression_method_from_path_extension(src));
        }

        if self.format.is_none() {
            self.format = detect_format_from_path_extension(src);
        }

        File::create(src)
            .map(BufWriter::new)
            .and_then(|writer| self.build_from_writer(writer))
    }

    /// Builds an alignment writer from a writer.
    ///
    /// If the format is not set, a default format is used. If the compression method is not set, a
    /// default one is determined by the format.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_util::alignment::{self, Format};
    /// let writer = alignment::writer::Builder::default().build_from_writer(io::sink());
    /// ```
    pub fn build_from_writer<W>(self, writer: W) -> io::Result<Writer>
    where
        W: Write + 'static,
    {
        let format = self.format.unwrap_or(Format::Sam);

        let compression_method = match self.compression_method {
            Some(compression_method) => compression_method,
            None => match format {
                Format::Sam | Format::Cram => None,
                Format::Bam => Some(CompressionMethod::Bgzf),
            },
        };

        let inner: Box<dyn sam::AlignmentWriter> = match (format, compression_method) {
            (Format::Sam, None) => Box::new(sam::Writer::new(writer)),
            (Format::Sam, Some(CompressionMethod::Bgzf)) => {
                Box::new(sam::Writer::new(bgzf::Writer::new(writer)))
            }
            (Format::Bam, None) => Box::new(bam::Writer::from(writer)),
            (Format::Bam, Some(CompressionMethod::Bgzf)) => Box::new(bam::Writer::new(writer)),
            (Format::Cram, None) => Box::new(
                cram::writer::Builder::default()
                    .set_reference_sequence_repository(self.reference_sequence_repository)
                    .set_block_content_encoder_map(self.block_content_encoder_map)
                    .build_with_writer(writer),
            ),
            (Format::Cram, Some(CompressionMethod::Bgzf)) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "CRAM cannot be bgzip-compressed",
                ));
            }
        };

        Ok(Writer { inner })
    }
}

fn detect_compression_method_from_path_extension<P>(path: P) -> Option<CompressionMethod>
where
    P: AsRef<Path>,
{
    match path.as_ref().extension().and_then(|ext| ext.to_str()) {
        Some("bam" | "gz" | "bgz") => Some(CompressionMethod::Bgzf),
        _ => None,
    }
}

fn detect_format_from_path_extension<P>(path: P) -> Option<Format>
where
    P: AsRef<Path>,
{
    let path = path.as_ref();

    match path.extension().and_then(|ext| ext.to_str()) {
        Some("sam") => Some(Format::Sam),
        Some("bam") => Some(Format::Bam),
        Some("cram") => Some(Format::Cram),
        Some("gz" | "bgz") => {
            let file_stem = path.file_stem().and_then(|stem| stem.to_str())?;
            file_stem.ends_with("sam").then_some(Format::Sam)
        }
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_compression_method_from_path_extension() {
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam"),
            None
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam.gz"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.sam.bgz"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.bam"),
            Some(CompressionMethod::Bgzf)
        );
        assert_eq!(
            detect_compression_method_from_path_extension("out.cram"),
            None
        );

        assert!(detect_format_from_path_extension("out.fa").is_none());
    }

    #[test]
    fn test_detect_format_from_path_extension() {
        assert_eq!(
            detect_format_from_path_extension("out.sam"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.sam.gz"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.sam.bgz"),
            Some(Format::Sam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.bam"),
            Some(Format::Bam)
        );
        assert_eq!(
            detect_format_from_path_extension("out.cram"),
            Some(Format::Cram)
        );

        assert!(detect_format_from_path_extension("out.fa").is_none());
    }
}
