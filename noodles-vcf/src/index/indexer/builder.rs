use noodles_csi::{self as csi, binning_index::index::header::ReferenceSequenceNames};

use super::{Indexer, Inner};
use crate::index::Format;

/// A VCF indexer builder.
#[derive(Default)]
pub struct Builder {
    format: Format,
}

impl Builder {
    /// Sets the VCF index format.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::{Format, Indexer};
    /// let builder = Indexer::builder().set_format(Format::Tabix);
    /// ```
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Builds a VCF indexer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_vcf::index::Indexer;
    /// let indexer = Indexer::builder().build();
    /// ```
    pub fn build(self) -> Indexer {
        let inner = match self.format {
            Format::Csi => Inner::Csi(csi::binning_index::Indexer::default()),
            Format::Tabix => Inner::Tabix(csi::binning_index::Indexer::default()),
        };

        Indexer {
            header: csi::binning_index::index::header::Builder::vcf().build(),
            reference_sequence_names: ReferenceSequenceNames::default(),
            inner,
        }
    }
}
