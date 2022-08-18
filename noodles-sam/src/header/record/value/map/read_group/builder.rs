use super::{Platform, ReadGroup};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference read group.
#[derive(Debug, Default)]
pub struct Builder {
    id: Option<String>,
    barcode: Option<String>,
    sequencing_center: Option<String>,
    description: Option<String>,
    produced_at: Option<String>,
    flow_order: Option<String>,
    key_sequence: Option<String>,
    library: Option<String>,
    program: Option<String>,
    predicted_median_insert_size: Option<i32>,
    platform: Option<Platform>,
    platform_model: Option<String>,
    platform_unit: Option<String>,
    sample: Option<String>,
}

impl map::Builder<ReadGroup> {
    /// Sets a read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    /// let read_group = Map::<ReadGroup>::builder().set_id("rg0").build()?;
    /// assert_eq!(read_group.id(), "rg0");
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_id<I>(mut self, id: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.id = Some(id.into());
        self
    }

    /// Sets a barcode sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_barcode("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.barcode(), Some("ACGT"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_barcode<I>(mut self, barcode: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.barcode = Some(barcode.into());
        self
    }

    /// Sets a sequencing center.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_sequencing_center("sc0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sequencing_center(), Some("sc0"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_sequencing_center<I>(mut self, sequencing_center: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.sequencing_center = Some(sequencing_center.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.description(), Some("noodles"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.description = Some(description.into());
        self
    }

    /// Sets a datetime of run.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_produced_at("2020-08-19T20:00:00Z")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.produced_at(), Some("2020-08-19T20:00:00Z"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_produced_at<I>(mut self, produced_at: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.produced_at = Some(produced_at.into());
        self
    }

    /// Sets a flow order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_flow_order("*")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.flow_order(), Some("*"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_flow_order<I>(mut self, flow_order: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.flow_order = Some(flow_order.into());
        self
    }

    /// Sets a key sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_key_sequence("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.key_sequence(), Some("ACGT"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_key_sequence<I>(mut self, key_sequence: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.key_sequence = Some(key_sequence.into());
        self
    }

    /// Sets a library.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_library("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.library(), Some("sample0"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_library<I>(mut self, library: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.library = Some(library.into());
        self
    }

    /// Sets a list of programs used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_program("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.program(), Some("noodles"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_program<I>(mut self, program: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.program = Some(program.into());
        self
    }

    /// Sets a predicted median insert size.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_predicted_median_insert_size(101)
    ///     .build()?;
    ///
    /// assert_eq!(read_group.predicted_median_insert_size(), Some(101));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_predicted_median_insert_size(mut self, predicted_median_insert_size: i32) -> Self {
        self.inner.predicted_median_insert_size = Some(predicted_median_insert_size);
        self
    }

    /// Sets a platform.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{
    ///     map::{read_group::Platform, ReadGroup},
    ///     Map,
    /// };
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_platform(Platform::Illumina)
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform(), Some(Platform::Illumina));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_platform(mut self, platform: Platform) -> Self {
        self.inner.platform = Some(platform);
        self
    }

    /// Sets a platform model.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_platform_model("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_model(), Some("noodles"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_platform_model<I>(mut self, platform_model: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.platform_model = Some(platform_model.into());
        self
    }

    /// Sets a platform unit.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_platform_unit("NDLS000.1")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_unit(), Some("NDLS000.1"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_platform_unit<I>(mut self, platform_unit: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.platform_unit = Some(platform_unit.into());
        self
    }

    /// Sets a sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_id("rg0")
    ///     .set_sample("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sample(), Some("sample0"));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_sample<I>(mut self, sample: I) -> Self
    where
        I: Into<String>,
    {
        self.inner.sample = Some(sample.into());
        self
    }
}

impl map::builder::Inner<ReadGroup> for Builder {
    fn build(self) -> Result<ReadGroup, BuildError> {
        let id = self.id.ok_or(BuildError::MissingField("ID"))?;

        Ok(ReadGroup {
            id,
            barcode: self.barcode,
            sequencing_center: self.sequencing_center,
            description: self.description,
            produced_at: self.produced_at,
            flow_order: self.flow_order,
            key_sequence: self.key_sequence,
            library: self.library,
            program: self.program,
            predicted_median_insert_size: self.predicted_median_insert_size,
            platform: self.platform,
            platform_model: self.platform_model,
            platform_unit: self.platform_unit,
            sample: self.sample,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let builder = Builder::default();

        assert!(builder.id.is_none());
        assert!(builder.barcode.is_none());
        assert!(builder.sequencing_center.is_none());
        assert!(builder.description.is_none());
        assert!(builder.produced_at.is_none());
        assert!(builder.flow_order.is_none());
        assert!(builder.key_sequence.is_none());
        assert!(builder.library.is_none());
        assert!(builder.program.is_none());
        assert!(builder.predicted_median_insert_size.is_none());
        assert!(builder.platform.is_none());
        assert!(builder.platform_model.is_none());
        assert!(builder.platform_unit.is_none());
        assert!(builder.sample.is_none());
    }
}
