use noodles_core::Region;
use serde::Deserialize;
use url::Url;

use super::{Class, Kind, Payload};
use crate::{Client, Error, Format, Response, response::Ticket};

/// A request builder.
pub struct Builder {
    client: Client,
    kind: Kind,
    id: String,
    payload: Payload,
}

impl Builder {
    pub fn new<I>(client: Client, kind: Kind, id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            client,
            kind,
            id: id.into(),
            payload: Payload::from(kind),
        }
    }

    pub fn set_format(mut self, format: Format) -> Self {
        *self.payload.format_mut() = format;
        self
    }

    pub fn set_class(mut self, class: Class) -> Self {
        *self.payload.class_mut() = Some(class);
        self
    }

    pub fn add_region(mut self, region: Region) -> Self {
        self.payload.regions_mut().push(region);
        self
    }

    pub async fn send(self) -> crate::Result<Response> {
        let endpoint = build_endpoint(self.client.base_url(), self.kind, &self.id)?;
        let mut request = self.client.http_client().post(endpoint);

        request = match self.payload.class() {
            Some(Class::Header) => {
                let mut payload = Payload::from(self.kind);
                *payload.class_mut() = Some(Class::Header);
                request.json(&payload)
            }
            None => request.json(&self.payload),
        };

        let response = request.send().await.map_err(Error::Request)?;

        if response.status().is_client_error() {
            let data: ErrorResponse = response.json().await.map_err(Error::Request)?;
            Err(Error::Response(data.htsget))
        } else {
            let data: TicketResponse = response.json().await.map_err(Error::Request)?;
            Ok(Response::new(self.client, self.id, data.htsget))
        }
    }
}

#[derive(Deserialize)]
pub struct TicketResponse {
    htsget: Ticket,
}

#[derive(Deserialize)]
pub struct ErrorResponse {
    htsget: crate::response::Error,
}

fn build_endpoint(base_url: &Url, kind: Kind, id: &str) -> crate::Result<Url> {
    let k = match kind {
        Kind::Reads => "reads",
        Kind::Variants => "variants",
    };

    base_url.join(&format!("{k}/{id}")).map_err(Error::Url)
}
