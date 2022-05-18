use noodles_core::Region;

use crate::{request, Response};

/// A reads endpoint builder.
pub struct Builder {
    inner: request::Builder,
}

impl Builder {
    pub(crate) fn new(inner: request::Builder) -> Self {
        Self { inner }
    }

    /// Adds a region to query.
    pub fn add_region(mut self, region: Region) -> Self {
        self.inner = self.inner.add_region(region);
        self
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Response> {
        self.inner.send().await
    }
}
