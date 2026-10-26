mod bounds;

use std::io;

use noodles_core::Position;

pub(crate) use self::bounds::Bounds;

const DEFAULT_BUF: &str = ".1.N....";
const MISSING: &str = ".";

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: String,
    pub(crate) bounds: Bounds,
}

impl Fields {
    pub(crate) fn reset(&mut self) {
        self.buf.clear();
        self.buf.push_str(DEFAULT_BUF);
        self.bounds = Bounds::default();
    }

    pub(super) fn reference_sequence_name(&self) -> &str {
        &self.buf[self.bounds.reference_sequence_name_range()]
    }

    pub(super) fn variant_start(&self) -> Option<io::Result<Position>> {
        const TELOMERE_START: &str = "0";

        match &self.buf[self.bounds.variant_start_range()] {
            TELOMERE_START => None,
            src => Some(
                src.parse::<usize>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    .and_then(|n| {
                        Position::try_from(n)
                            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
                    }),
            ),
        }
    }

    pub(super) fn ids(&self) -> &str {
        match &self.buf[self.bounds.ids_range()] {
            MISSING => "",
            src => src,
        }
    }

    pub(super) fn reference_bases(&self) -> &str {
        &self.buf[self.bounds.reference_bases_range()]
    }

    pub(super) fn alternate_bases(&self) -> &str {
        match &self.buf[self.bounds.alternate_bases_range()] {
            MISSING => "",
            src => src,
        }
    }

    pub(super) fn quality_score(&self) -> Option<&str> {
        match &self.buf[self.bounds.quality_score_range()] {
            MISSING => None,
            src => Some(src),
        }
    }

    pub(super) fn filters(&self) -> &str {
        match &self.buf[self.bounds.filters_range()] {
            MISSING => "",
            src => src,
        }
    }

    pub(super) fn info(&self) -> &str {
        match &self.buf[self.bounds.info_range()] {
            MISSING => "",
            src => src,
        }
    }

    pub(super) fn samples_keys(&self) -> &str {
        &self.buf[self.bounds.samples_keys_range()]
    }

    pub(super) fn samples_values(&self) -> &str {
        &self.buf[self.bounds.samples_values_range()]
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            buf: String::from(DEFAULT_BUF),
            bounds: Bounds::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reset() {
        let mut fields = Fields {
            buf: String::from(".1.N....GT0|0"),
            bounds: Bounds {
                samples_keys_end: 10,
                ..Default::default()
            },
        };

        fields.reset();

        assert_eq!(fields, Fields::default());
    }

    #[test]
    fn test_samples() {
        let fields = Fields::default();
        assert!(fields.samples_keys().is_empty());
        assert!(fields.samples_values().is_empty());

        let fields = Fields {
            buf: String::from(".1.N....GT0|0"),
            bounds: Bounds {
                samples_keys_end: 10,
                ..Default::default()
            },
        };
        assert_eq!(fields.samples_keys(), "GT");
        assert_eq!(fields.samples_values(), "0|0");

        let fields = Fields {
            buf: String::from(".1.N......"),
            bounds: Bounds {
                samples_keys_end: 9,
                ..Default::default()
            },
        };
        assert_eq!(fields.samples_keys(), ".");
        assert_eq!(fields.samples_values(), ".");
    }
}
