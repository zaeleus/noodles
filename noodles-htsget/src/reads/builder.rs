use serde::Deserialize;

use super::Reads;
use crate::{Client, Error, Ticket};

/// A reads endpoint builder.
pub struct Builder {
    client: Client,
    id: String,
    reference_name: Option<String>,
    start: Option<u32>,
    end: Option<u32>,
}

impl Builder {
    pub(crate) fn new<I>(client: Client, id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            client,
            id: id.into(),
            reference_name: None,
            start: None,
            end: None,
        }
    }

    /// Sets the reference name.
    pub fn set_reference_name<N>(mut self, reference_name: N) -> Self
    where
        N: Into<String>,
    {
        self.reference_name = Some(reference_name.into());
        self
    }

    /// Sets the start position.
    ///
    /// This is 0-based, inclusive.
    pub fn set_start(mut self, start: u32) -> Self {
        self.start = Some(start);
        self
    }

    /// Sets the end position.
    ///
    /// This is 0-based, exclusive.
    pub fn set_end(mut self, end: u32) -> Self {
        self.end = Some(end);
        self
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Reads> {
        let endpoint = self
            .client
            .base_url()
            .join(&format!("reads/{}", self.id))
            .map_err(Error::Url)?;
        let mut request = self.client.http_client().get(endpoint);

        let mut query = Vec::new();

        if let Some(reference_name) = self.reference_name {
            query.push(("referenceName", reference_name));
        }

        if let Some(start) = self.start {
            query.push(("start", start.to_string()));
        }

        if let Some(end) = self.end {
            query.push(("end", end.to_string()));
        }

        request = request.query(&query);

        let response = request.send().await.map_err(Error::Request)?;
        let data: TicketResponse = response.json().await.map_err(Error::Request)?;

        Ok(Reads::new(self.client, self.id, data.htsget))
    }
}

#[derive(Deserialize)]
pub struct TicketResponse {
    htsget: Ticket,
}
