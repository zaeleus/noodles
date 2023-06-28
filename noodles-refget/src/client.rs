use url::Url;

use super::sequence;

/// A refget client.
#[derive(Clone, Debug)]
pub struct Client {
    http_client: reqwest::Client,
    base_url: Url,
}

impl Client {
    /// Creates a new refget client.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_refget as refget;
    /// let client = refget::Client::new("https://localhost/".parse()?);
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn new(base_url: Url) -> Self {
        Self {
            http_client: reqwest::Client::new(),
            base_url,
        }
    }

    /// Creates a refget client with the given HTTP client.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_refget as refget;
    /// let http_client = reqwest::Client::new();
    /// let base_url = "https://localhost/".parse()?;
    /// let client = refget::Client::with_http_client(http_client, base_url);
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

    /// Creates a sequence request for the given ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_refget as refget;
    /// let client = refget::Client::new("https://localhost/".parse()?);
    /// let sequence_builder = client.sequence("d7eba311421bbc9d3ada44709dd61534");
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn sequence<I>(&self, id: I) -> sequence::Builder
    where
        I: Into<String>,
    {
        sequence::Builder::new(self.clone(), id)
    }

    /// Creates a sequence metadata request for the given ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_refget as refget;
    /// let client = refget::Client::new("https://localhost/".parse()?);
    /// let sequence_builder = client.sequence_metadata("d7eba311421bbc9d3ada44709dd61534");
    /// # Ok::<_, url::ParseError>(())
    /// ```
    pub fn sequence_metadata<I>(&self, id: I) -> sequence::metadata::Builder
    where
        I: Into<String>,
    {
        sequence::metadata::Builder::new(self.clone(), id)
    }
}
