use noodles_core::Region;
use serde::Deserialize;

use crate::{Client, Error, Response, Ticket};

/// A variants endpoint builder.
pub struct Builder {
    client: Client,
    id: String,
    region: Option<Region>,
}

impl Builder {
    pub(crate) fn new<I>(client: Client, id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            client,
            id: id.into(),
            region: None,
        }
    }

    /// Sets the region to query.
    pub fn set_region(mut self, region: Region) -> Self {
        self.region = Some(region);
        self
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Response> {
        use crate::resolve_interval;

        let endpoint = self
            .client
            .base_url()
            .join(&format!("variants/{}", self.id))
            .map_err(Error::Url)?;
        let mut request = self.client.http_client().get(endpoint);

        if let Some(region) = self.region {
            let mut query = Vec::with_capacity(3);

            query.push(("referenceName", region.name().into()));

            let (resolved_start, resolved_end) = resolve_interval(region.interval());

            if let Some(position) = resolved_start {
                let start = u32::try_from(position).map_err(|_| Error::Input)?;
                query.push(("start", start.to_string()));
            }

            if let Some(position) = resolved_end {
                let end = u32::try_from(position).map_err(|_| Error::Input)?;
                query.push(("end", end.to_string()));
            }

            request = request.query(&query);
        }

        let response = request.send().await.map_err(Error::Request)?;
        let data: TicketResponse = response.json().await.map_err(Error::Request)?;

        Ok(Response::new(self.client, self.id, data.htsget))
    }
}

#[derive(Deserialize)]
pub struct TicketResponse {
    htsget: Ticket,
}
