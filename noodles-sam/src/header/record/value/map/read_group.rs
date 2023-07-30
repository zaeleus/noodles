//! SAM header record read group map value.

mod builder;
pub mod platform;
pub(crate) mod tag;

pub use self::platform::Platform;
pub(crate) use self::tag::Tag;

use std::fmt;

use self::{builder::Builder, tag::StandardTag};
use super::{Inner, Map};

/// A SAM header record read group map value.
///
/// A read group typically defines the set of reads that came from the same run on a sequencing
/// instrument.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ReadGroup {
    pub(crate) barcode: Option<String>,
    pub(crate) sequencing_center: Option<String>,
    pub(crate) description: Option<String>,
    pub(crate) produced_at: Option<String>,
    pub(crate) flow_order: Option<String>,
    pub(crate) key_sequence: Option<String>,
    pub(crate) library: Option<String>,
    pub(crate) program: Option<String>,
    pub(crate) predicted_median_insert_size: Option<i32>,
    pub(crate) platform: Option<Platform>,
    pub(crate) platform_model: Option<String>,
    pub(crate) platform_unit: Option<String>,
    pub(crate) sample: Option<String>,
}

impl Inner for ReadGroup {
    type StandardTag = StandardTag;
    type Builder = Builder;
}

impl Map<ReadGroup> {
    /// Returns the barcode sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.barcode().is_none());
    /// ```
    pub fn barcode(&self) -> Option<&str> {
        self.inner.barcode.as_deref()
    }

    /// Returns the sequencing center.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.sequencing_center().is_none());
    /// ```
    pub fn sequencing_center(&self) -> Option<&str> {
        self.inner.sequencing_center.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.description().is_none());
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.inner.description.as_deref()
    }

    /// Returns the datatime of run.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.produced_at().is_none());
    /// ```
    pub fn produced_at(&self) -> Option<&str> {
        self.inner.produced_at.as_deref()
    }

    /// Returns the flow order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.flow_order().is_none());
    /// ```
    pub fn flow_order(&self) -> Option<&str> {
        self.inner.flow_order.as_deref()
    }

    /// Returns the key sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.key_sequence().is_none());
    /// ```
    pub fn key_sequence(&self) -> Option<&str> {
        self.inner.key_sequence.as_deref()
    }

    /// Returns the library.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.library().is_none());
    /// ```
    pub fn library(&self) -> Option<&str> {
        self.inner.library.as_deref()
    }

    /// Returns the programs used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.program().is_none());
    /// ```
    pub fn program(&self) -> Option<&str> {
        self.inner.program.as_deref()
    }

    /// Returns the predicted median insert size, rounded to the nearest integer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.predicted_median_insert_size().is_none());
    /// ```
    pub fn predicted_median_insert_size(&self) -> Option<i32> {
        self.inner.predicted_median_insert_size
    }

    /// Returns the platform used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.platform().is_none());
    /// ```
    pub fn platform(&self) -> Option<Platform> {
        self.inner.platform
    }

    /// Returns the platform model.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.platform_model().is_none());
    /// ```
    pub fn platform_model(&self) -> Option<&str> {
        self.inner.platform_model.as_deref()
    }

    /// Returns the platform unit.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.platform_unit().is_none());
    /// ```
    pub fn platform_unit(&self) -> Option<&str> {
        self.inner.platform_unit.as_deref()
    }

    /// Returns the sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::default();
    /// assert!(read_group.sample().is_none());
    /// ```
    pub fn sample(&self) -> Option<&str> {
        self.inner.sample.as_deref()
    }
}

impl fmt::Display for Map<ReadGroup> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(barcode) = self.barcode() {
            write!(f, "\t{}:{barcode}", tag::BARCODE)?;
        }

        if let Some(sequencing_center) = self.sequencing_center() {
            write!(f, "\t{}:{sequencing_center}", tag::SEQUENCING_CENTER)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{description}", tag::DESCRIPTION)?;
        }

        if let Some(produced_at) = self.produced_at() {
            write!(f, "\t{}:{produced_at}", tag::PRODUCED_AT)?;
        }

        if let Some(flow_order) = self.flow_order() {
            write!(f, "\t{}:{flow_order}", tag::FLOW_ORDER)?;
        }

        if let Some(key_sequence) = self.key_sequence() {
            write!(f, "\t{}:{key_sequence}", tag::KEY_SEQUENCE)?;
        }

        if let Some(library) = self.library() {
            write!(f, "\t{}:{library}", tag::LIBRARY)?;
        }

        if let Some(program) = self.program() {
            write!(f, "\t{}:{program}", tag::PROGRAM)?;
        }

        if let Some(predicted_median_insert_size) = self.predicted_median_insert_size() {
            write!(
                f,
                "\t{}:{predicted_median_insert_size}",
                tag::PREDICTED_MEDIAN_INSERT_SIZE
            )?;
        }

        if let Some(platform) = self.platform() {
            write!(f, "\t{}:{platform}", tag::PLATFORM)?;
        }

        if let Some(platform_model) = self.platform_model() {
            write!(f, "\t{}:{platform_model}", tag::PLATFORM_MODEL)?;
        }

        if let Some(platform_unit) = self.platform_unit() {
            write!(f, "\t{}:{platform_unit}", tag::PLATFORM_UNIT)?;
        }

        if let Some(sample) = self.sample() {
            write!(f, "\t{}:{sample}", tag::SAMPLE)?;
        }

        super::fmt_display_other_fields(f, self.other_fields())?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::record::value::map::builder::BuildError;

    #[test]
    fn test_fmt() -> Result<(), BuildError> {
        let read_group = Map::<ReadGroup>::builder()
            .set_program("noodles")
            .set_platform(Platform::Illumina)
            .build()?;

        assert_eq!(read_group.to_string(), "\tPG:noodles\tPL:ILLUMINA");

        Ok(())
    }
}
