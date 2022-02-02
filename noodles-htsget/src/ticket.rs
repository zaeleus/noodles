use std::collections::HashMap;

use serde::Deserialize;
use url::Url;

use super::Format;

#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub(crate) struct BlockUrl {
    url: Url,
    #[serde(default)]
    headers: HashMap<String, String>,
    class: Option<String>,
}

impl BlockUrl {
    pub fn url(&self) -> &Url {
        &self.url
    }

    pub fn headers(&self) -> &HashMap<String, String> {
        &self.headers
    }
}

#[derive(Clone, Debug, Deserialize, Eq, PartialEq)]
pub(crate) struct Ticket {
    format: Format,
    urls: Vec<BlockUrl>,
    md5: Option<String>,
}

impl Ticket {
    pub fn urls(&self) -> &[BlockUrl] {
        &self.urls
    }
}
