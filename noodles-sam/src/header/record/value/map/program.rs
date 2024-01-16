//! SAM header record program map value.

mod builder;
pub mod tag;

pub(crate) use self::tag::Tag;

use self::builder::Builder;
use super::Inner;

// A SAM header record program map value.
///
/// A program describes any program that created, viewed, or mutated a SAM file.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Program;

impl Inner for Program {
    type StandardTag = tag::Standard;
    type Builder = Builder;
}
