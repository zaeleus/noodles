//! SAM header read group and fields.

mod builder;
pub mod platform;
pub mod tag;

pub use self::{builder::Builder, platform::Platform, tag::Tag};

use std::{collections::HashMap, convert::TryFrom, error, fmt, num};

use super::{
    record::{self, value::Fields},
    Record,
};

/// A SAM header read group.
///
/// A read group typically defines the set of reads that came from the same run on a sequencing
/// instrument. The read group ID is guaranteed to be set.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReadGroup {
    id: String,
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

impl ReadGroup {
    /// Creates a SAM header read group builder.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let builder = ReadGroup::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Creates a read group with an ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn new<I>(id: I) -> Self
    where
        I: Into<String>,
    {
        Self {
            id: id.into(),
            barcode: None,
            sequencing_center: None,
            description: None,
            produced_at: None,
            flow_order: None,
            key_sequence: None,
            library: None,
            program: None,
            predicted_median_insert_size: None,
            platform: None,
            platform_model: None,
            platform_unit: None,
            sample: None,
            fields: HashMap::new(),
        }
    }

    /// Returns the read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Returns a mutable reference to the read group ID.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    ///
    /// let mut read_group = ReadGroup::new("rg0");
    /// assert_eq!(read_group.id(), "rg0");
    ///
    /// *read_group.id_mut() = String::from("rg1");
    /// assert_eq!(read_group.id(), "rg1");
    /// ```
    pub fn id_mut(&mut self) -> &mut String {
        &mut self.id
    }

    /// Returns the barcode sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.barcode().is_none());
    /// ```
    pub fn barcode(&self) -> Option<&str> {
        self.barcode.as_deref()
    }

    /// Returns the sequencing center.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.sequencing_center().is_none());
    /// ```
    pub fn sequencing_center(&self) -> Option<&str> {
        self.sequencing_center.as_deref()
    }

    /// Returns the description.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.description().is_none());
    /// ```
    pub fn description(&self) -> Option<&str> {
        self.description.as_deref()
    }

    /// Returns the datatime of run.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.produced_at().is_none());
    /// ```
    pub fn produced_at(&self) -> Option<&str> {
        self.produced_at.as_deref()
    }

    /// Returns the flow order.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.flow_order().is_none());
    /// ```
    pub fn flow_order(&self) -> Option<&str> {
        self.flow_order.as_deref()
    }

    /// Returns the key sequence.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.key_sequence().is_none());
    /// ```
    pub fn key_sequence(&self) -> Option<&str> {
        self.key_sequence.as_deref()
    }

    /// Returns the library.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.library().is_none());
    /// ```
    pub fn library(&self) -> Option<&str> {
        self.library.as_deref()
    }

    /// Returns the programs used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.program().is_none());
    /// ```
    pub fn program(&self) -> Option<&str> {
        self.program.as_deref()
    }

    /// Returns the predicted median insert size.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.predicted_median_insert_size().is_none());
    /// ```
    pub fn predicted_median_insert_size(&self) -> Option<i32> {
        self.predicted_median_insert_size
    }

    /// Returns the platform used.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.platform().is_none());
    /// ```
    pub fn platform(&self) -> Option<Platform> {
        self.platform
    }

    /// Returns the platform model.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.platform_model().is_none());
    /// ```
    pub fn platform_model(&self) -> Option<&str> {
        self.platform_model.as_deref()
    }

    /// Returns the platform unit.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.platform_unit().is_none());
    /// ```
    pub fn platform_unit(&self) -> Option<&str> {
        self.platform_unit.as_deref()
    }

    /// Returns the sample.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::ReadGroup;
    /// let read_group = ReadGroup::new("rg0");
    /// assert!(read_group.sample().is_none());
    /// ```
    pub fn sample(&self) -> Option<&str> {
        self.sample.as_deref()
    }

    /// Returns the raw fields of the read group.
    ///
    /// This includes any field that is not specially handled by the structure itself. For example,
    /// this will not include the ID field, as it is parsed and available as [`Self::id`].
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::header::{read_group::Tag, ReadGroup};
    ///
    /// let read_group = ReadGroup::builder()
    ///     .set_id("rg0")
    ///     .insert(Tag::Other(String::from("zn")), String::from("noodles"))
    ///     .build();
    ///
    /// let fields = read_group.fields();
    /// assert_eq!(fields.len(), 1);
    /// assert_eq!(
    ///     fields.get(&Tag::Other(String::from("zn"))),
    ///     Some(&String::from("noodles"))
    /// );
    ///
    /// assert_eq!(fields.get(&Tag::Id), None);
    /// assert_eq!(read_group.id(), "rg0");
    /// ```
    pub fn fields(&self) -> &HashMap<Tag, String> {
        &self.fields
    }
}

impl fmt::Display for ReadGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", record::Kind::ReadGroup)?;
        write!(f, "\t{}:{}", Tag::Id, self.id)?;

        if let Some(barcode) = self.barcode() {
            write!(f, "\t{}:{}", Tag::Barcode, barcode)?;
        }

        if let Some(sequencing_center) = self.sequencing_center() {
            write!(f, "\t{}:{}", Tag::SequencingCenter, sequencing_center)?;
        }

        if let Some(description) = self.description() {
            write!(f, "\t{}:{}", Tag::Description, description)?;
        }

        if let Some(produced_at) = self.produced_at() {
            write!(f, "\t{}:{}", Tag::ProducedAt, produced_at)?;
        }

        if let Some(flow_order) = self.flow_order() {
            write!(f, "\t{}:{}", Tag::FlowOrder, flow_order)?;
        }

        if let Some(key_sequence) = self.key_sequence() {
            write!(f, "\t{}:{}", Tag::KeySequence, key_sequence)?;
        }

        if let Some(library) = self.library() {
            write!(f, "\t{}:{}", Tag::Library, library)?;
        }

        if let Some(program) = self.program() {
            write!(f, "\t{}:{}", Tag::Program, program)?;
        }

        if let Some(predicted_median_insert_size) = self.predicted_median_insert_size() {
            write!(
                f,
                "\t{}:{}",
                Tag::PredictedMedianInsertSize,
                predicted_median_insert_size
            )?;
        }

        if let Some(platform) = self.platform() {
            write!(f, "\t{}:{}", Tag::Platform, platform)?;
        }

        if let Some(platform_model) = self.platform_model() {
            write!(f, "\t{}:{}", Tag::PlatformModel, platform_model)?;
        }

        if let Some(platform_unit) = self.platform_unit() {
            write!(f, "\t{}:{}", Tag::PlatformUnit, platform_unit)?;
        }

        if let Some(sample) = self.sample() {
            write!(f, "\t{}:{}", Tag::Sample, sample)?;
        }

        for (tag, value) in &self.fields {
            write!(f, "\t{}:{}", tag, value)?;
        }

        Ok(())
    }
}

/// An error returned when a raw SAM header read group fails to parse.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TryFromRecordError {
    /// The record is invalid.
    InvalidRecord,
    /// A required tag is missing.
    MissingRequiredTag(Tag),
    /// A tag is invalid.
    InvalidTag(tag::ParseError),
    /// The predicted median insert size is invalid.
    InvalidPredictedMedianInsertSize(num::ParseIntError),
    /// The platform is invalid.
    InvalidPlatform(platform::ParseError),
}

impl error::Error for TryFromRecordError {}

impl fmt::Display for TryFromRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRecord => f.write_str("invalid record"),
            Self::MissingRequiredTag(tag) => write!(f, "missing required tag: {:?}", tag),
            Self::InvalidTag(e) => write!(f, "invalid tag: {}", e),
            Self::InvalidPredictedMedianInsertSize(e) => {
                write!(f, "invalid predicted median insert size: {}", e)
            }
            Self::InvalidPlatform(e) => write!(f, "invalid platform: {}", e),
        }
    }
}

impl TryFrom<Record> for ReadGroup {
    type Error = TryFromRecordError;

    fn try_from(record: Record) -> Result<Self, Self::Error> {
        match record.into() {
            (record::Kind::ReadGroup, record::Value::Map(fields)) => parse_map(fields),
            _ => Err(TryFromRecordError::InvalidRecord),
        }
    }
}

fn parse_map(raw_fields: Fields) -> Result<ReadGroup, TryFromRecordError> {
    let mut builder = ReadGroup::builder();
    let mut id = None;

    for (raw_tag, value) in raw_fields {
        let tag = raw_tag.parse().map_err(TryFromRecordError::InvalidTag)?;

        builder = match tag {
            Tag::Id => {
                id = Some(value);
                builder
            }
            Tag::Barcode => builder.set_barcode(value),
            Tag::SequencingCenter => builder.set_sequencing_center(value),
            Tag::Description => builder.set_description(value),
            Tag::ProducedAt => builder.set_produced_at(value),
            Tag::FlowOrder => builder.set_flow_order(value),
            Tag::KeySequence => builder.set_key_sequence(value),
            Tag::Library => builder.set_library(value),
            Tag::Program => builder.set_program(value),
            Tag::PredictedMedianInsertSize => {
                let predicted_median_insert_size = value
                    .parse()
                    .map_err(TryFromRecordError::InvalidPredictedMedianInsertSize)?;

                builder.set_predicted_median_insert_size(predicted_median_insert_size)
            }
            Tag::Platform => {
                let platform = value.parse().map_err(TryFromRecordError::InvalidPlatform)?;
                builder.set_platform(platform)
            }
            Tag::PlatformModel => builder.set_platform_model(value),
            Tag::PlatformUnit => builder.set_platform_unit(value),
            Tag::Sample => builder.set_sample(value),
            Tag::Other(..) => builder.insert(tag, value),
        }
    }

    if let Some(i) = id {
        builder = builder.set_id(i);
    } else {
        return Err(TryFromRecordError::MissingRequiredTag(Tag::Id));
    }

    Ok(builder.build())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fmt() {
        let read_group = ReadGroup::builder()
            .set_id("rg0")
            .set_program("noodles")
            .build();

        assert_eq!(read_group.to_string(), "@RG\tID:rg0\tPG:noodles");
    }

    #[test]
    fn test_try_from_record_for_read_group_with_invalid_record() {
        let record = Record::new(
            record::Kind::Comment,
            record::Value::String(String::from("noodles-sam")),
        );

        assert_eq!(
            ReadGroup::try_from(record),
            Err(TryFromRecordError::InvalidRecord)
        );
    }

    #[test]
    fn test_try_from_record_for_read_group_with_no_id(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReadGroup,
            record::Value::try_from_iter(vec![("PG", "noodles")])?,
        );

        assert_eq!(
            ReadGroup::try_from(record),
            Err(TryFromRecordError::MissingRequiredTag(Tag::Id))
        );

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_read_group_with_invalid_predicted_median_insert_size(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReadGroup,
            record::Value::try_from_iter(vec![("ID", "rg0"), ("PI", "unknown")])?,
        );

        assert!(matches!(
            ReadGroup::try_from(record),
            Err(TryFromRecordError::InvalidPredictedMedianInsertSize(_))
        ));

        Ok(())
    }

    #[test]
    fn test_try_from_record_for_read_group_with_invalid_platform(
    ) -> Result<(), record::value::TryFromIteratorError> {
        let record = Record::new(
            record::Kind::ReadGroup,
            record::Value::try_from_iter(vec![("ID", "rg0"), ("PL", "unknown")])?,
        );

        assert!(matches!(
            ReadGroup::try_from(record),
            Err(TryFromRecordError::InvalidPlatform(_))
        ));

        Ok(())
    }
}
