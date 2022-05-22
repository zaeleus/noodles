use std::ops::RangeBounds;

use noodles_core::{region::Interval, Position};

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

        let response = request.send().await.map_err(Error::Request)?;
        let sequence = response.bytes().await.map_err(Error::Request)?;

        Ok(Sequence::new(self.client, self.id, sequence))
    }
}

// Resolves a 1-based interval as a 0-based [start, end) interval.
fn resolve_interval<B>(interval: B) -> (Option<usize>, Option<usize>)
where
    B: RangeBounds<Position>,
{
    use std::ops::Bound;

    let start = match interval.start_bound() {
        Bound::Included(position) => Some(usize::from(*position) - 1),
        Bound::Excluded(position) => Some(usize::from(*position)),
        Bound::Unbounded => None,
    };

    let end = match interval.end_bound() {
        Bound::Included(position) => Some(usize::from(*position)),
        Bound::Excluded(position) => Some(usize::from(*position) - 1),
        Bound::Unbounded => None,
    };

    (start, end)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resolve_interval() -> std::result::Result<(), noodles_core::position::TryFromIntError> {
        let start = Position::try_from(8)?;
        let end = Position::try_from(13)?;

        assert_eq!(resolve_interval(start..end), (Some(7), Some(12)));
        assert_eq!(resolve_interval(start..), (Some(7), None));
        assert_eq!(resolve_interval(..), (None, None));
        assert_eq!(resolve_interval(start..=end), (Some(7), Some(13)));
        assert_eq!(resolve_interval(..end), (None, Some(12)));
        assert_eq!(resolve_interval(..=end), (None, Some(13)));

        Ok(())
    }
}
