use url::Url;

use super::{reads, request, request::Kind, variants};

/// A htsget client.
#[derive(Clone, Debug)]
pub struct Client {
    http_client: reqwest::Client,
    base_url: Url,
}

impl Client {
    /// Creates an htsget client with a default HTTP client.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_htsget as htsget;
    /// let base_url = "https://localhost/".parse()?;
    /// let client = htsget::Client::new(base_url);
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn new(base_url: Url) -> Self {
        Self::with_http_client(reqwest::Client::new(), base_url)
    }

    /// Creates a htsget client with the given HTTP client.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_htsget as htsget;
    /// let http_client = reqwest::Client::new();
    /// let base_url = "https://localhost/".parse()?;
    /// let client = htsget::Client::with_http_client(http_client, base_url);
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn with_http_client(http_client: reqwest::Client, base_url: Url) -> Self {
        Self {
            http_client,
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
        let builder = request::Builder::new(self.clone(), Kind::Reads, id);
        reads::Builder::new(builder)
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
        let builder = request::Builder::new(self.clone(), Kind::Variants, id);
        variants::Builder::new(builder)
    }
}
