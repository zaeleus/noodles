use url::Url;

use crate::variants;

use super::reads;

/// A htsget client.
#[derive(Clone, Debug)]
pub struct Client {
    http_client: reqwest::Client,
    base_url: Url,
}

impl Client {
    /// Creates a new htsget client.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_htsget as htsget;
    /// let client = htsget::Client::new("https://localhost/".parse()?);
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn new(base_url: Url) -> Self {
        Self {
            http_client: reqwest::Client::new(),
            base_url,
        }
    }

    pub(crate) fn http_client(&self) -> &reqwest::Client {
        &self.http_client
    }

    pub(crate) fn base_url(&self) -> &Url {
        &self.base_url
    }

    /// Creates a reads request for the given ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_htsget as htsget;
    /// let client = htsget::Client::new("https://localhost/".parse()?);
    /// let reads = client.reads("NDLS0001");
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn reads<I>(&self, id: I) -> reads::Builder
    where
        I: Into<String>,
    {
        reads::Builder::new(self.clone(), id)
    }

    /// Creates a variants request for the given ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_htsget as htsget;
    /// let client = htsget::Client::new("https://localhost/".parse()?);
    /// let variants = client.variants("NDLS0001");
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn variants<I>(&self, id: I) -> variants::Builder
    where
        I: Into<String>,
    {
        variants::Builder::new(self.clone(), id)
    }
}
