//! SAM header read group builder.

use std::{collections::HashMap, error, fmt};

use super::{Platform, ReadGroup, Tag};

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
    fields: HashMap<Tag, String>,
}

/// An error returned when a SAM header read group fails to build.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BuildError {
    /// The ID is missing.
    MissingId,
}

impl error::Error for BuildError {}

impl fmt::Display for BuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingId => f.write_str("missing ID"),
        }
    }
}

impl Builder {
    /// Sets a read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::builder().set_id("rg0").build()?;
    /// assert_eq!(read_group.id(), "rg0");
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_id<I>(mut self, id: I) -> Self
    where
        I: Into<String>,
    {
        self.id = Some(id.into());
        self
    }

    /// Sets a barcode sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_barcode("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.barcode(), Some("ACGT"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_barcode<I>(mut self, barcode: I) -> Self
    where
        I: Into<String>,
    {
        self.barcode = Some(barcode.into());
        self
    }

    /// Sets a sequencing center.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_sequencing_center("sc0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sequencing_center(), Some("sc0"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_sequencing_center<I>(mut self, sequencing_center: I) -> Self
    where
        I: Into<String>,
    {
        self.sequencing_center = Some(sequencing_center.into());
        self
    }

    /// Sets a description.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_description("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.description(), Some("noodles"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_description<I>(mut self, description: I) -> Self
    where
        I: Into<String>,
    {
        self.description = Some(description.into());
        self
    }

    /// Sets a datetime of run.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_produced_at("2020-08-19T20:00:00Z")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.produced_at(), Some("2020-08-19T20:00:00Z"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_produced_at<I>(mut self, produced_at: I) -> Self
    where
        I: Into<String>,
    {
        self.produced_at = Some(produced_at.into());
        self
    }

    /// Sets a flow order.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_flow_order("*")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.flow_order(), Some("*"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_flow_order<I>(mut self, flow_order: I) -> Self
    where
        I: Into<String>,
    {
        self.flow_order = Some(flow_order.into());
        self
    }

    /// Sets a key sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_key_sequence("ACGT")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.key_sequence(), Some("ACGT"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_key_sequence<I>(mut self, key_sequence: I) -> Self
    where
        I: Into<String>,
    {
        self.key_sequence = Some(key_sequence.into());
        self
    }

    /// Sets a library.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_library("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.library(), Some("sample0"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_library<I>(mut self, library: I) -> Self
    where
        I: Into<String>,
    {
        self.library = Some(library.into());
        self
    }

    /// Sets a list of programs used.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_program("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.program(), Some("noodles"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_program<I>(mut self, program: I) -> Self
    where
        I: Into<String>,
    {
        self.program = Some(program.into());
        self
    }

    /// Sets a predicted median insert size.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_predicted_median_insert_size(101)
    ///     .build()?;
    ///
    /// assert_eq!(read_group.predicted_median_insert_size(), Some(101));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_predicted_median_insert_size(mut self, predicted_median_insert_size: i32) -> Self {
        self.predicted_median_insert_size = Some(predicted_median_insert_size);
        self
    }

    /// Sets a platform.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::{read_group::Platform, ReadGroup};
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_platform(Platform::Illumina)
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform(), Some(Platform::Illumina));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_platform(mut self, platform: Platform) -> Self {
        self.platform = Some(platform);
        self
    }

    /// Sets a platform model.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_platform_model("noodles")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_model(), Some("noodles"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_platform_model<I>(mut self, platform_model: I) -> Self
    where
        I: Into<String>,
    {
        self.platform_model = Some(platform_model.into());
        self
    }

    /// Sets a platform unit.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_platform_unit("NDLS000.1")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.platform_unit(), Some("NDLS000.1"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_platform_unit<I>(mut self, platform_unit: I) -> Self
    where
        I: Into<String>,
    {
        self.platform_unit = Some(platform_unit.into());
        self
    }

    /// Sets a sample.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .set_sample("sample0")
    ///     .build()?;
    ///
    /// assert_eq!(read_group.sample(), Some("sample0"));
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn set_sample<I>(mut self, sample: I) -> Self
    where
        I: Into<String>,
    {
        self.sample = Some(sample.into());
        self
    }

    /// Inserts a tag-raw value pair.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::{read_group::Tag, ReadGroup};
    ///
    /// let zn = Tag::Other(String::from("zn"));
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .insert(zn.clone(), String::from("noodles"))
    ///     .build()?;
    ///
    /// assert_eq!(
    ///     read_group.fields().get(&zn),
    ///     Some(&String::from("noodles"))
    /// );
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn insert<I>(mut self, tag: Tag, value: I) -> Self
    where
        I: Into<String>,
    {
        self.fields.insert(tag, value.into());
        self
    }

    /// Builds a reference sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// # use noodles_sam::header::read_group::builder;
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::builder().set_id("rg0").build()?;
    /// # Ok::<(), builder::BuildError>(())
    /// ```
    pub fn build(self) -> Result<ReadGroup, BuildError> {
        let id = self.id.ok_or(BuildError::MissingId)?;

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
            fields: self.fields,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build() {
        assert_eq!(Builder::default().build(), Err(BuildError::MissingId));
    }
}
