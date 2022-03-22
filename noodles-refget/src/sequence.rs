//! Sequence endpoint.

mod builder;
pub mod metadata;

use bytes::Bytes;
use serde::Deserialize;

use crate::{Client, Error};

pub use self::{builder::Builder, metadata::Metadata};

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

    /// Returns metadata related to the sequence.
    pub async fn metadata(&self) -> crate::Result<Metadata> {
        let endpoint = self
            .client
            .base_url()
            .join(&format!("sequence/{}/metadata", self.id))
            .map_err(Error::Url)?;

        let response = self
            .client
            .http_client()
            .get(endpoint)
            .send()
            .await
            .map_err(Error::Request)?;

        response
            .json()
            .await
            .map(|data: MetadataResponse| data.metadata)
            .map_err(Error::Request)
    }
}

#[derive(Deserialize)]
struct MetadataResponse {
    metadata: Metadata,
}
