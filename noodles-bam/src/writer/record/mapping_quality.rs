use bytes::BufMut;
use noodles_sam::{self as sam, alignment::record::MappingQuality};

pub fn put_mapping_quality<B>(dst: &mut B, mapping_quality: Option<MappingQuality>)
where
    B: BufMut,
{
    use sam::alignment::record::mapping_quality::MISSING;
    let mapq = mapping_quality.map(u8::from).unwrap_or(MISSING);
    dst.put_u8(mapq);
}
