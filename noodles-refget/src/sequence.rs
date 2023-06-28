//! Sequence endpoint.

mod builder;
pub mod metadata;
pub mod service;

use bytes::Bytes;

use crate::Client;

pub use self::{builder::Builder, metadata::Metadata, service::Service};

/// A response from the sequence endpoint.
#[derive(Debug)]
pub struct Sequence {
    client: Client,
    id: String,
    sequence: Bytes,
}

impl Sequence {
    pub(crate) fn new(client: Client, id: String, sequence: Bytes) -> Self {
        Self {
            client,
            id,
            sequence,
        }
    }

    /// Returns the sequence.
    pub fn sequence(&self) -> Bytes {
        self.sequence.clone()
    }

    /// Builds a request to get the metadata related to the sequence.
    pub fn metadata(&self) -> metadata::Builder {
        metadata::Builder::new(self.client.clone(), &self.id)
    }
}
