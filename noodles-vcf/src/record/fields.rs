mod bounds;

use std::io;

use noodles_core::Position;

pub(crate) use self::bounds::Bounds;

const MISSING: &str = ".";

#[derive(Clone, Eq, PartialEq)]
pub(crate) struct Fields {
    pub(crate) buf: String,
    pub(crate) bounds: Bounds,
}

impl Fields {
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

    pub(super) fn samples(&self) -> &str {
        &self.buf[self.bounds.samples_range()]
    }
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            buf: String::from(".1.N...."),
            bounds: Bounds::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_samples() {
        let fields = Fields::default();
        assert!(fields.samples().is_empty());

        let fields = Fields {
            buf: String::from(".1.N....GT\t0|0"),
            ..Default::default()
        };
        assert_eq!(fields.samples(), "GT\t0|0");

        let fields = Fields {
            buf: String::from(".1.N.....\t."),
            ..Default::default()
        };
        assert_eq!(fields.samples(), ".\t.");
    }
}
