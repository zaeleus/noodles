//! Reads endpoint.

mod builder;

pub use self::builder::Builder;

use bytes::Bytes;
use futures::Stream;

use super::{Client, Ticket};

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
        use super::chunks::chunks;
        chunks(&self.client, self.ticket.urls())
    }
}
