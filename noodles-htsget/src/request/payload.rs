mod regions;

use noodles_core::Region;
use serde::Serialize;

use self::regions::Regions;
use super::{Class, Kind};
use crate::Format;

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct Payload {
    format: Format,

    #[serde(skip_serializing_if = "Option::is_none")]
    class: Option<Class>,

    #[serde(skip_serializing_if = "Regions::is_empty")]
    regions: Regions,
}

impl Payload {
    pub fn format_mut(&mut self) -> &mut Format {
        &mut self.format
    }

    pub fn class(&self) -> Option<Class> {
        self.class
    }

    pub fn class_mut(&mut self) -> &mut Option<Class> {
        &mut self.class
    }

    pub fn regions_mut(&mut self) -> &mut Vec<Region> {
        &mut self.regions.0
    }
}

impl From<Kind> for Payload {
    fn from(kind: Kind) -> Self {
        let format = match kind {
            Kind::Reads => Format::Bam,
            Kind::Variants => Format::Vcf,
        };

        Self {
            format,
            class: None,
            regions: Regions::default(),
        }
    }
}

#[cfg(test)]
mod tests {
    use serde_test::{Token, assert_ser_tokens};

    use super::*;

    #[test]
    fn test_serialize() {
        let payload = Payload::from(Kind::Reads);

        assert_ser_tokens(
            &payload,
            &[
                Token::Struct {
                    name: "Payload",
                    len: 1,
                },
                Token::Str("format"),
                Token::UnitVariant {
                    name: "Format",
                    variant: "BAM",
                },
                Token::StructEnd,
            ],
        );

        let mut payload = Payload::from(Kind::Reads);
        *payload.class_mut() = Some(Class::Header);

        assert_ser_tokens(
            &payload,
            &[
                Token::Struct {
                    name: "Payload",
                    len: 2,
                },
                Token::Str("format"),
                Token::UnitVariant {
                    name: "Format",
                    variant: "BAM",
                },
                Token::Str("class"),
                Token::Some,
                Token::UnitVariant {
                    name: "Class",
                    variant: "header",
                },
                Token::StructEnd,
            ],
        );

        let mut payload = Payload::from(Kind::Reads);
        payload.regions_mut().push(Region::new("sq0", ..));

        assert_ser_tokens(
            &payload,
            &[
                Token::Struct {
                    name: "Payload",
                    len: 2,
                },
                Token::Str("format"),
                Token::UnitVariant {
                    name: "Format",
                    variant: "BAM",
                },
                Token::Str("regions"),
                Token::Seq { len: Some(1) },
                Token::Map { len: None },
                Token::Str("referenceName"),
                Token::Str("sq0"),
                Token::MapEnd,
                Token::SeqEnd,
                Token::StructEnd,
            ],
        );
    }
}
