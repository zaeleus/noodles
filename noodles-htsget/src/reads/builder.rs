use noodles_core::Region;

use super::Format;
use crate::{
    Response,
    request::{self, Class},
};

/// A reads endpoint builder.
pub struct Builder {
    inner: request::Builder,
}

impl Builder {
    pub(crate) fn new(inner: request::Builder) -> Self {
        Self { inner }
    }

    /// Sets the data format.
    pub fn set_format(mut self, format: Format) -> Self {
        self.inner = self.inner.set_format(format.into());
        self
    }

    /// Sets the data stream filter.
    ///
    /// Setting this discards all other options upon send.
    pub fn set_class(mut self, class: Class) -> Self {
        self.inner = self.inner.set_class(class);
        self
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
