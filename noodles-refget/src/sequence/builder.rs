use noodles_core::region::Interval;

use crate::{Client, Error, Sequence};

/// A sequence endpoint builder.
pub struct Builder {
    client: Client,
    id: String,
    interval: Option<Interval>,
}

impl Builder {
    pub(crate) fn new<I>(client: Client, id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            client,
            id: id.into(),
            interval: None,
        }
    }

    /// Sets the interval to query.
    pub fn set_interval<I>(mut self, interval: I) -> Self
    where
        I: Into<Interval>,
    {
        self.interval = Some(interval.into());
        self
    }

    /// Sends the request.
    pub async fn send(self) -> crate::Result<Sequence> {
        let endpoint = self
            .client
            .base_url()
            .join(&format!("sequence/{}", self.id))
            .map_err(Error::Url)?;

        let mut request = self.client.http_client().get(endpoint);

        if let Some(interval) = self.interval {
            let mut query = Vec::new();

            let (resolved_start, resolved_end) = resolve_interval(interval);

            if let Some(start) = resolved_start {
                query.push(("start", start.to_string()));
            }

            if let Some(end) = resolved_end {
                query.push(("end", end.to_string()));
            }

            request = request.query(&query);
        }

        let response = request.send().await?.error_for_status()?;
        let sequence = response.bytes().await?;

        Ok(Sequence::new(self.client, self.id, sequence))
    }
}

fn resolve_interval<I>(interval: I) -> (Option<usize>, Option<usize>)
where
    I: Into<Interval>,
{
    let interval = interval.into();
    let start = interval.start().map(|position| usize::from(position) - 1);
    let end = interval.end().map(usize::from);
    (start, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_interval() -> std::result::Result<(), noodles_core::position::TryFromIntError> {
        use noodles_core::Position;

        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        assert_eq!(resolve_interval(start..=end), (Some(7), Some(13)));
        assert_eq!(resolve_interval(start..), (Some(7), None));
        assert_eq!(resolve_interval(..=end), (None, Some(13)));
        assert_eq!(resolve_interval(..), (None, None));

        Ok(())
    }
}
