use std::io;

use super::{Data, RecordBuf};
use crate::{alignment::Record, Header};

impl RecordBuf {
    /// Converts an alignment record to a buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_sam::{self as sam, alignment::RecordBuf};
    ///
    /// let header = sam::Header::default();
    /// let record = sam::Record::default();
    ///
    /// let record_buf = RecordBuf::try_from_alignment_record(&header, &record)?;
    ///
    /// assert_eq!(record_buf, RecordBuf::default());
    /// # Ok::<_, std::io::Error>(())
    /// ```
    pub fn try_from_alignment_record<R>(header: &Header, record: &R) -> io::Result<Self>
    where
        R: Record,
    {
        let mut record_buf = RecordBuf::default();

        *record_buf.name_mut() = record.name().map(|name| name.into());
        *record_buf.flags_mut() = record.flags()?;
        *record_buf.reference_sequence_id_mut() =
            record.reference_sequence_id(header).transpose()?;
        *record_buf.alignment_start_mut() = record.alignment_start().transpose()?;
        *record_buf.mapping_quality_mut() = record.mapping_quality().transpose()?;
        *record_buf.cigar_mut() = record.cigar().iter().collect::<io::Result<_>>()?;
        *record_buf.mate_reference_sequence_id_mut() =
            record.mate_reference_sequence_id(header).transpose()?;
        *record_buf.mate_alignment_start_mut() = record.mate_alignment_start().transpose()?;
        *record_buf.template_length_mut() = record.template_length()?;
        *record_buf.sequence_mut() = record.sequence().iter().collect::<Vec<_>>().into();
        *record_buf.quality_scores_mut() =
            record.quality_scores().iter().collect::<Vec<_>>().into();

        let mut data_buf = Data::default();

        for result in record.data().iter() {
            let (tag, value) = result?;
            data_buf.insert(tag, value.try_into()?);
        }

        *record_buf.data_mut() = data_buf;

        Ok(record_buf)
    }
}
