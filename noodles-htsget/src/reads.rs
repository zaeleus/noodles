//! Reads endpoint.

mod builder;

pub use self::builder::Builder;

use std::pin::Pin;

use bytes::Bytes;
use futures::{stream, Stream, TryStreamExt};

use super::{ticket::BlockUrl, Client, Error, Ticket};

/// A response from the reads endpoint.
#[derive(Debug)]
pub struct Reads {
    client: Client,
    id: String,
    ticket: Ticket,
}

impl Reads {
    pub(crate) fn new(client: Client, id: String, ticket: Ticket) -> Self {
        Self { client, id, ticket }
    }

    /// Returns the record ID.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the data from the ticket URLs.
    pub fn chunks(&self) -> impl Stream<Item = crate::Result<Bytes>> + '_ {
        Box::pin(
            stream::try_unfold(
                (&self.client, self.ticket.urls(), 0),
                |(client, urls, i)| async move {
                    match urls.get(i) {
                        Some(url) => {
                            let st = resolve_data(client, url).await;
                            Ok(Some((st, (client, urls, i + 1))))
                        }
                        None => Ok(None),
                    }
                },
            )
            .try_flatten(),
        )
    }
}

async fn resolve_data(
    client: &Client,
    block_url: &BlockUrl,
) -> Pin<Box<dyn Stream<Item = crate::Result<Bytes>>>> {
    const PREFIX: &str = "data:application/octet-stream;base64,";

    let url = block_url.url();

    if url.scheme() == "data" {
        if let Some(encoded_data) = url.as_str().strip_prefix(PREFIX) {
            match base64::decode(encoded_data) {
                Ok(data) => Box::pin(stream::once(async { Ok(Bytes::from(data)) })),
                Err(e) => Box::pin(stream::once(async { Err(Error::Decode(e)) })),
            }
        } else {
            Box::pin(stream::once(async { Err(Error::InvalidDataUrl) }))
        }
    } else {
        let mut request = client.http_client().get(url.clone());

        for (key, value) in block_url.headers() {
            request = request.header(key, value);
        }

        match request.send().await.map_err(Error::Request) {
            Ok(response) => Box::pin(response.bytes_stream().map_err(Error::Request)),
            Err(e) => Box::pin(stream::once(async { Err(e) })),
        }
    }
}
