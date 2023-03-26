//! htsget response.

mod error;
pub(crate) mod ticket;

pub use self::error::Error;
pub(crate) use self::ticket::Ticket;

use bytes::Bytes;
use futures::Stream;

use super::Client;

/// An htsget response.
#[derive(Debug)]
pub struct Response {
    client: Client,
    id: String,
    ticket: Ticket,
}

impl Response {
    pub(crate) fn new(client: Client, id: String, ticket: Ticket) -> Self {
        Self { client, id, ticket }
    }

    /// Returns the record ID.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns the data from the ticket URLs.
    pub fn chunks(&self) -> impl Stream<Item = crate::Result<Bytes>> + '_ {
        use super::chunks::chunks;
        chunks(&self.client, self.ticket.urls())
    }
}
