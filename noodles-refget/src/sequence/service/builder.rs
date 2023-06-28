use serde::Deserialize;

use super::Service;
use crate::{Client, Error};

/// A sequence metadata endpoint builder.
pub struct Builder {
    client: Client,
}

impl Builder {
    pub(crate) fn new(client: Client) -> Self {
        Self { client }
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Service> {
        let endpoint = self
            .client
            .base_url()
            .join("sequence/service-info")
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
            .map(|data: ServiceInfoResponse| data.serivce)
            .map_err(Error::Request)
    }
}

#[derive(Deserialize)]
struct ServiceInfoResponse {
    serivce: Service,
}
