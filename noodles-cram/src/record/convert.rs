use std::{io, str};

use noodles_core::Position;
use noodles_sam as sam;

use super::{Features, Flags, QualityScores, Record, Sequence};

impl Record {
    /// Converts an alignment record to a CRAM record.
    pub fn try_from_alignment_record<R>(header: &sam::Header, record: &R) -> io::Result<Self>
    where
        R: sam::alignment::Record,
    {
        let mut builder = Self::builder();

        let bam_flags = record.flags().try_to_u16().map(sam::record::Flags::from)?;
        builder = builder.set_bam_flags(bam_flags);

        let mut flags = Flags::default();

        if let Some(reference_sequence_id) = record.reference_sequence_id(header) {
            let id = reference_sequence_id.try_to_usize()?;
            builder = builder.set_reference_sequence_id(id);
        }

        builder = builder.set_read_length(record.sequence().len());

        if let Some(alignment_start) = record.alignment_start() {
            let position = alignment_start.try_to_usize().and_then(|n| {
                Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

            builder = builder.set_alignment_start(position);
        }

        let mut data = alignment_record_data_to_data_buf(record.data())?;

        if let Some(read_group_id) = get_read_group_id(header.read_groups(), &data)? {
            builder = builder.set_read_group_id(read_group_id);
        }

        if let Some(name) = record.name() {
            let name = sam::alignment::record_buf::Name::from(name.as_bytes());
            builder = builder.set_name(name);
        }

        // next mate bit flags

        if let Some(mate_reference_sequence_id) = record.mate_reference_sequence_id(header) {
            let id = mate_reference_sequence_id.try_to_usize()?;
            builder = builder.set_next_fragment_reference_sequence_id(id);
        }

        if let Some(mate_alignment_start) = record.mate_alignment_start() {
            let position = mate_alignment_start.try_to_usize().and_then(|n| {
                Position::try_from(n).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

            builder = builder.set_next_mate_alignment_start(position);
        }

        let template_length = record.template_length().try_to_i32()?;
        builder = builder.set_template_size(template_length);

        // distance to next fragment

        if !data.is_empty() {
            use sam::record::data::field::tag;
            data.remove(&tag::READ_GROUP);
            builder = builder.set_tags(data);
        }

        let raw_bases: Vec<_> = record.sequence().iter().collect();
        let bases = Sequence::from(raw_bases);

        let quality_scores = if record.quality_scores().is_empty() {
            QualityScores::default()
        } else {
            if bam_flags.is_unmapped() {
                flags.insert(Flags::QUALITY_SCORES_STORED_AS_ARRAY);
            }

            let scores: Vec<_> = record.quality_scores().iter().collect();
            QualityScores::from(scores)
        };

        if !bam_flags.is_unmapped() {
            let cigar = alignment_record_cigar_to_cigar_buf(record.cigar(header))?;
            let features = Features::from_cigar(flags, &cigar, &bases, &quality_scores);
            builder = builder.set_features(features);
        }

        if let Some(mapping_quality) = record.mapping_quality() {
            let mapping_quality = mapping_quality.try_to_u8().and_then(|n| {
                sam::record::MappingQuality::try_from(n)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
            })?;

            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !record.quality_scores().is_empty() {
            if bam_flags.is_unmapped() {
                flags.insert(Flags::QUALITY_SCORES_STORED_AS_ARRAY);
            }

            let scores: Vec<_> = record.quality_scores().iter().collect();
            let quality_scores = QualityScores::from(scores);
            builder = builder.set_quality_scores(quality_scores);
        }

        Ok(builder
            .set_flags(flags)
            .set_bases(bases)
            .set_quality_scores(quality_scores)
            .build())
    }

    /// Converts this CRAM record to an alignment record.
    pub fn try_into_alignment_record(
        self,
        header: &sam::Header,
    ) -> io::Result<sam::alignment::RecordBuf> {
        let mut builder = sam::alignment::RecordBuf::builder();

        if let Some(read_name) = self.name {
            builder = builder.set_name(read_name);
        }

        builder = builder.set_flags(self.bam_bit_flags);

        if let Some(reference_sequence_id) = self.reference_sequence_id {
            builder = builder.set_reference_sequence_id(reference_sequence_id);
        }

        if let Some(alignment_start) = self.alignment_start {
            builder = builder.set_alignment_start(alignment_start);
        }

        if let Some(mapping_quality) = self.mapping_quality {
            builder = builder.set_mapping_quality(mapping_quality);
        }

        if !self.bam_bit_flags.is_unmapped() {
            let cigar = self
                .features
                .try_into_cigar(self.read_length)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

            builder = builder.set_cigar(cigar);
        }

        if let Some(mate_reference_sequence_id) = self.next_fragment_reference_sequence_id {
            builder = builder.set_mate_reference_sequence_id(mate_reference_sequence_id);
        }

        if let Some(mate_alignment_start) = self.next_mate_alignment_start {
            builder = builder.set_mate_alignment_start(mate_alignment_start);
        }

        builder = builder
            .set_template_length(self.template_size)
            .set_sequence(self.bases)
            .set_quality_scores(self.quality_scores);

        let mut data = self.tags;
        maybe_insert_read_group(&mut data, header.read_groups(), self.read_group_id)?;
        builder = builder.set_data(data);

        Ok(builder.build())
    }
}

fn alignment_record_cigar_to_cigar_buf<C>(cigar: C) -> io::Result<sam::record::Cigar>
where
    C: sam::alignment::record::Cigar,
{
    use sam::record::cigar::Op;

    cigar
        .iter()
        .map(|result| result.map(|(kind, len)| Op::new(kind, len)))
        .collect()
}

fn alignment_record_data_to_data_buf<D>(data: D) -> io::Result<sam::record::Data>
where
    D: sam::alignment::record::Data,
{
    use sam::{
        alignment::record::data::field::{value::Array, Value},
        record::data::field::{
            value::Array as ArrayBuf,
            value::{Character, Hex},
            Tag, Value as ValueBuf,
        },
    };

    fn value_to_value_buf(value: Value<'_>) -> io::Result<ValueBuf> {
        match value {
            Value::Character(c) => Character::try_from(c)
                .map(ValueBuf::Character)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
            Value::Int8(n) => Ok(ValueBuf::Int8(n)),
            Value::UInt8(n) => Ok(ValueBuf::UInt8(n)),
            Value::Int16(n) => Ok(ValueBuf::Int16(n)),
            Value::UInt16(n) => Ok(ValueBuf::UInt16(n)),
            Value::Int32(n) => Ok(ValueBuf::Int32(n)),
            Value::UInt32(n) => Ok(ValueBuf::UInt32(n)),
            Value::Float(n) => Ok(ValueBuf::Float(n)),
            Value::String(s) => str::from_utf8(s)
                .map(|t| ValueBuf::String(t.into()))
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
            Value::Hex(s) => Hex::try_from(s)
                .map(ValueBuf::Hex)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
            Value::Array(Array::Int8(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::Int8(vs))),
            Value::Array(Array::UInt8(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::UInt8(vs))),
            Value::Array(Array::Int16(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::Int16(vs))),
            Value::Array(Array::UInt16(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::UInt16(vs))),
            Value::Array(Array::Int32(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::Int32(vs))),
            Value::Array(Array::UInt32(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::UInt32(vs))),
            Value::Array(Array::Float(values)) => values
                .iter()
                .collect::<Result<_, _>>()
                .map(|vs| ValueBuf::Array(ArrayBuf::Float(vs))),
        }
    }

    let mut buf = sam::record::Data::default();

    for result in data.iter() {
        let (raw_tag, raw_value) = result?;

        let tag =
            Tag::try_from(raw_tag).map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

        let value = value_to_value_buf(raw_value)?;

        buf.insert(tag, value);
    }

    Ok(buf)
}

fn get_read_group_id(
    read_groups: &sam::header::ReadGroups,
    data: &sam::record::Data,
) -> io::Result<Option<usize>> {
    use sam::record::data::field::tag;

    let Some(rg_value) = data.get(&tag::READ_GROUP) else {
        return Ok(None);
    };

    let read_group_name = rg_value.as_str().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "invalid read group field value",
        )
    })?;

    read_groups
        .get_index_of(read_group_name)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read group name"))
        .map(Some)
}

fn maybe_insert_read_group(
    data: &mut sam::record::Data,
    read_groups: &sam::header::ReadGroups,
    read_group_id: Option<usize>,
) -> io::Result<()> {
    use sam::record::data::field::{tag, Value};

    if let Some(id) = read_group_id {
        let name = read_groups
            .get_index(id)
            .map(|(name, _)| name)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "invalid read group ID"))?;

        data.insert(tag::READ_GROUP, Value::String(name.into()));
    }

    Ok(())
}
