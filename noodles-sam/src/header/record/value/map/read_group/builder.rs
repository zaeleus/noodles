use super::{Platform, ReadGroup};
use crate::header::record::value::map::{self, builder::BuildError};

/// A SAM header reference read group.
#[derive(Debug, Default)]
pub struct Builder {
    barcode: Option<Vec<u8>>,
    sequencing_center: Option<Vec<u8>>,
    description: Option<Vec<u8>>,
    produced_at: Option<Vec<u8>>,
    flow_order: Option<Vec<u8>>,
    key_sequence: Option<Vec<u8>>,
    library: Option<Vec<u8>>,
    program: Option<Vec<u8>>,
    predicted_median_insert_size: Option<i32>,
    platform: Option<Platform>,
    platform_model: Option<Vec<u8>>,
    platform_unit: Option<Vec<u8>>,
    sample: Option<Vec<u8>>,
}

impl map::Builder<ReadGroup> {
    /// Sets a barcode sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::record::value::{map::ReadGroup, Map};
    ///
    /// let read_group = Map::<ReadGroup>::builder()
    ///     .set_barcode("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.barcode(), Some(&b"ACGT"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_barcode<I>(mut self, barcode: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_sequencing_center("sc0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sequencing_center(), Some(&b"sc0"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_sequencing_center<I>(mut self, sequencing_center: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.description(), Some(&b"noodles"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_produced_at("2020-08-19T20:00:00Z")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.produced_at(), Some(&b"2020-08-19T20:00:00Z"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_produced_at<I>(mut self, produced_at: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_flow_order("*")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.flow_order(), Some(&b"*"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_flow_order<I>(mut self, flow_order: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_key_sequence("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.key_sequence(), Some(&b"ACGT"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_key_sequence<I>(mut self, key_sequence: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_library("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.library(), Some(&b"sample0"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_library<I>(mut self, library: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_program("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.program(), Some(&b"noodles"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_program<I>(mut self, program: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_platform_model("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_model(), Some(&b"noodles"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_platform_model<I>(mut self, platform_model: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_platform_unit("NDLS000.1")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_unit(), Some(&b"NDLS000.1"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_platform_unit<I>(mut self, platform_unit: I) -> Self
    where
        I: Into<Vec<u8>>,
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
    ///     .set_sample("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sample(), Some(&b"sample0"[..]));
    /// # Ok::<_, noodles_sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn set_sample<I>(mut self, sample: I) -> Self
    where
        I: Into<Vec<u8>>,
    {
        self.inner.sample = Some(sample.into());
        self
    }
}

impl map::builder::Inner<ReadGroup> for Builder {
    fn build(self) -> Result<ReadGroup, BuildError> {
        Ok(ReadGroup {
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
