use noodles_csi as csi;

use super::{Indexer, Inner};
use crate::index::Format;

/// A BAM indexer builder.
#[derive(Default)]
pub struct Builder {
    format: Format,
}

impl Builder {
    /// Sets the BAM index format.
    pub fn set_format(mut self, format: Format) -> Self {
        self.format = format;
        self
    }

    /// Builds a BAM indexer.
    pub fn build(self) -> Indexer {
        let inner = match self.format {
            Format::Bai => Inner::Bai(csi::binning_index::Indexer::default()),
            Format::Csi => Inner::Csi(csi::binning_index::Indexer::default()),
        };

        Indexer { inner }
    }
}
