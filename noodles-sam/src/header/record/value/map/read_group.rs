//! SAM header record read group map value.

mod builder;
pub mod platform;
pub mod tag;

pub use self::platform::Platform;
pub(crate) use self::tag::Tag;

use self::builder::Builder;
use super::Inner;

/// A SAM header record read group map value.
///
/// A read group typically defines the set of reads that came from the same run on a sequencing
/// instrument.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadGroup;

impl Inner for ReadGroup {
    type StandardTag = tag::Standard;
    type Builder = Builder;
}
