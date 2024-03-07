mod key;
mod values;

use std::io::{self, Write};

use noodles_vcf::{self as vcf, header::StringMaps, variant::record::Samples};

use self::{
    key::write_key,
    values::{write_genotype_values, write_values},
};

pub fn write_samples<W, S>(
    writer: &mut W,
    header: &vcf::Header,
    string_maps: &StringMaps,
    samples: S,
) -> io::Result<()>
where
    W: Write,
    S: Samples,
{
    use noodles_vcf::variant::record_buf::samples::keys::key;

    for (i, result) in samples.column_names(header).enumerate() {
        let key = result?;

        write_key(writer, string_maps.strings(), key)?;

        let rows: Vec<_> = samples.iter().collect();
        let mut values = Vec::new();

        for sample in rows.iter() {
            let value = sample.get_index(header, i).transpose()?.unwrap_or_default();
            values.push(value);
        }

        let format = header.formats().get(key).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "missing FORMAT header record")
        })?;

        if key == key::GENOTYPE {
            write_genotype_values(writer, &values)?;
        } else {
            write_values(writer, format, &values)?;
        }

        drop(values);
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use noodles_vcf::header::StringMaps;

    use super::*;

    #[test]
    fn test_write_samples() -> Result<(), Box<dyn std::error::Error>> {
        use vcf::{
            header::record::value::Map,
            variant::record_buf::samples::{keys::key, sample::Value},
        };

        let header = vcf::Header::builder()
            .add_format(
                key::CONDITIONAL_GENOTYPE_QUALITY,
                Map::from(key::CONDITIONAL_GENOTYPE_QUALITY),
            )
            .add_format(key::READ_DEPTH, Map::from(key::READ_DEPTH))
            .add_sample_name("sample0")
            .add_sample_name("sample1")
            .build();

        let string_maps = StringMaps::try_from(&header)?;

        let genotypes = vcf::variant::record_buf::Samples::new(
            vcf::variant::record_buf::samples::Keys::try_from(vec![
                String::from(key::CONDITIONAL_GENOTYPE_QUALITY),
                String::from(key::READ_DEPTH),
            ])?,
            vec![
                vec![Some(Value::from(13)), Some(Value::from(5))],
                vec![Some(Value::from(8))],
            ],
        );

        let mut buf = Vec::new();
        write_samples(&mut buf, &header, &string_maps, &genotypes)?;

        let expected = [
            0x11, // string string map index type = Some(Type::Int(1))
            0x01, // string string map index = Some(1) (GQ)
            0x11, // GQ value type = Some(Type::Int8(1))
            0x0d, // GQ[0] = Some(13)
            0x08, // GQ[1] = Some(8)
            0x11, // string string map index type = Some(Type::Int(1))
            0x02, // string string map index = Some(2) (DP)
            0x11, // DP value type = Some(Type::Int8(1))
            0x05, // DP[0] = Some(5)
            0x80, // DP[1] = None
        ];

        assert_eq!(buf, expected);

        Ok(())
    }
}
