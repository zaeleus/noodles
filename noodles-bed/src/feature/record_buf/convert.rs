use std::io;

use super::{other_fields::Value, OtherFields, RecordBuf};
use crate::feature::Record;

impl RecordBuf<3> {
    /// Converts a BED3+ feature record to a BED3+ record buffer.
    pub fn try_from_feature_record<R>(record: &R) -> io::Result<Self>
    where
        R: Record<3>,
    {
        let mut builder = RecordBuf::<3>::builder()
            .set_reference_sequence_name(record.reference_sequence_name())
            .set_feature_start(record.feature_start()?);

        if let Some(feature_end) = record.feature_end().transpose()? {
            builder = builder.set_feature_end(feature_end);
        }

        let values: Vec<_> = record.other_fields().iter().map(Value::from).collect();
        let other_fields = OtherFields::from(values);
        builder = builder.set_other_fields(other_fields);

        Ok(builder.build())
    }
}

impl RecordBuf<4> {
    /// Converts a BED4+ feature record to a BED4+ record buffer.
    pub fn try_from_feature_record<R>(record: &R) -> io::Result<Self>
    where
        R: Record<4>,
    {
        let mut builder = RecordBuf::<4>::builder()
            .set_reference_sequence_name(record.reference_sequence_name())
            .set_feature_start(record.feature_start()?);

        if let Some(feature_end) = record.feature_end().transpose()? {
            builder = builder.set_feature_end(feature_end);
        }

        if let Some(Some(name)) = record.name() {
            builder = builder.set_name(name);
        }

        let values: Vec<_> = record.other_fields().iter().map(Value::from).collect();
        let other_fields = OtherFields::from(values);
        builder = builder.set_other_fields(other_fields);

        Ok(builder.build())
    }
}

impl RecordBuf<5> {
    /// Converts a BED4+ feature record to a BED4+ record buffer.
    pub fn try_from_feature_record<R>(record: &R) -> io::Result<Self>
    where
        R: Record<5>,
    {
        let mut builder = RecordBuf::<5>::builder()
            .set_reference_sequence_name(record.reference_sequence_name())
            .set_feature_start(record.feature_start()?);

        if let Some(feature_end) = record.feature_end().transpose()? {
            builder = builder.set_feature_end(feature_end);
        }

        if let Some(Some(name)) = record.name() {
            builder = builder.set_name(name);
        }

        if let Some(score) = record.score().transpose()? {
            builder = builder.set_score(score);
        }

        let values: Vec<_> = record.other_fields().iter().map(Value::from).collect();
        let other_fields = OtherFields::from(values);
        builder = builder.set_other_fields(other_fields);

        Ok(builder.build())
    }
}
