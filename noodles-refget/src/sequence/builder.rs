use crate::{Client, Error, Sequence};

/// A sequence endpoint builder.
pub struct Builder {
    client: Client,
    id: String,
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
            start: None,
            end: None,
        }
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
    pub async fn send(self) -> crate::Result<Sequence> {
        let endpoint = self
            .client
            .base_url()
            .join(&format!("sequence/{}", self.id))
            .map_err(Error::Url)?;

        let mut request = self.client.http_client().get(endpoint);

        let mut query = Vec::new();

        if let Some(start) = self.start {
            query.push(("start", start.to_string()));
        }

        if let Some(end) = self.end {
            query.push(("end", end.to_string()));
        }

        request = request.query(&query);

        let response = request.send().await.map_err(Error::Request)?;
        let sequence = response.bytes().await.map_err(Error::Request)?;

        Ok(Sequence::new(self.client, self.id, sequence))
    }
}
