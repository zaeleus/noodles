use super::Metadata;
use crate::{Client, Error};

use serde::Deserialize;

/// A sequence metadata endpoint builder.
pub struct Builder {
    client: Client,
    id: String,
}

impl Builder {
    pub(crate) fn new<I>(client: Client, id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            client,
            id: id.into(),
        }
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Metadata> {
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
            .await?
            .error_for_status()?;

        let metadata = response
            .json()
            .await
            .map(|data: MetadataResponse| data.metadata)?;

        Ok(metadata)
    }
}

#[derive(Deserialize)]
struct MetadataResponse {
    metadata: Metadata,
}
