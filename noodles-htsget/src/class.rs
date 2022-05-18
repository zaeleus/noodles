use serde::Deserialize;

/// The class of data.
#[derive(Clone, Copy, Debug, Deserialize, Eq, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Class {
    Header,
}
