use std::io::{self, Write};

use noodles_vcf::{self as vcf, variant::record::Filters as _};

use crate::header::string_maps::StringStringMap;

pub(super) fn write_filters<W>(
    writer: &mut W,
    string_string_map: &StringStringMap,
    filters: &vcf::variant::record_buf::Filters,
) -> io::Result<()>
where
    W: Write,
{
    use crate::record::codec::encoder::string_map::write_string_map_indices;

    let indices: Vec<_> = filters
        .iter()
        .map(|id| {
            string_string_map.get_index_of(id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("filter missing from string map: {id}"),
                )
            })
        })
        .collect::<Result<_, _>>()?;

    write_string_map_indices(writer, &indices)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::header::StringMaps;

    #[test]
    fn test_write_filters() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::{
            header::record::value::{map::Filter, Map},
            variant::record_buf::Filters,
        };

        fn t(
            buf: &mut Vec<u8>,
            string_map: &StringStringMap,
            filters: &Filters,
            expected: &[u8],
        ) -> io::Result<()> {
            buf.clear();
            write_filters(buf, string_map, filters)?;
            assert_eq!(buf, expected);
            Ok(())
        }

        let header = vcf::Header::builder()
            .add_filter("PASS", Map::<Filter>::pass())
            .add_filter(
                "s50",
                Map::<Filter>::new("Less than 50% of samples have data"),
            )
            .add_filter("q10", Map::<Filter>::new("Quality below 10"))
            .build();

        let string_maps = StringMaps::try_from(&header)?;

        let mut buf = Vec::new();

        let filters = Filters::default();
        t(&mut buf, string_maps.strings(), &filters, &[0x00])?;

        let filters = Filters::pass();
        t(&mut buf, string_maps.strings(), &filters, &[0x11, 0x00])?;

        let filters = [String::from("q10")].into_iter().collect();
        t(&mut buf, string_maps.strings(), &filters, &[0x11, 0x02])?;

        let filters = [String::from("q10"), String::from("s50")]
            .into_iter()
            .collect();
        t(
            &mut buf,
            string_maps.strings(),
            &filters,
            &[0x21, 0x02, 0x01],
        )?;

        Ok(())
    }
}
