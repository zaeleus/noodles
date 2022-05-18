use noodles_core::Region;

use crate::{request, Response};

/// A variants endpoint builder.
pub struct Builder {
    inner: request::Builder,
}

impl Builder {
    pub(crate) fn new(inner: request::Builder) -> Self {
        Self { inner }
    }

    /// Sets the region to query.
    pub fn set_region(mut self, region: Region) -> Self {
        self.inner = self.inner.set_region(region);
        self
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Response> {
        self.inner.send().await
    }
}
