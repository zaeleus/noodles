mod field;

use std::io;

use noodles_vcf::{self as vcf, variant::record::info::field::Value};

use self::field::read_field;
use crate::header::string_maps::StringStringMap;

/// BCF record info.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Info {
    buf: Vec<u8>,
    field_count: usize,
}

impl Info {
    /// Converts BCF record info to VCF record info.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf as bcf;
    /// use noodles_vcf as vcf;
    ///
    /// let bcf_info = bcf::record::Info::default();
    /// let header = vcf::Header::default();
    /// let string_maps = bcf::header::StringMaps::default();
    ///
    /// let vcf_info = bcf_info.try_into_vcf_record_info(&header, string_maps.strings())?;
    /// assert!(vcf_info.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_info<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_string_map: &'h StringStringMap,
    ) -> io::Result<vcf::record::Info> {
        let mut info = vcf::record::Info::default();

        for result in self.iter(header, string_string_map) {
            let (key, value) = result?;
            let value = value.map(|v| v.try_into()).transpose()?;
            info.insert(key.into(), value);
        }

        Ok(info)
    }

    /// Creates an info map by wrapping the given buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Info;
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    /// ```
    pub fn new(buf: Vec<u8>, field_count: usize) -> Self {
        Self { buf, field_count }
    }

    /// Returns the number of info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Info;
    /// let info = Info::default();
    /// assert_eq!(info.len(), 0);
    /// ```
    pub fn len(&self) -> usize {
        self.field_count
    }

    /// Returns whether there are any info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Info;
    /// let info = Info::default();
    /// assert!(info.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Removes all fields from the info map.
    ///
    /// This does not affect the capacity of the map.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::record::Info;
    /// let mut info = Info::default();
    /// info.clear();
    /// assert!(info.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.buf.clear();
        self.set_field_count(0);
    }

    /// Returns the value with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::{header::StringMaps, record::Info};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map, Map},
    ///     record::info::field::key,
    ///     variant::record::info::field::Value,
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(key::ALLELE_COUNT, Map::<map::Info>::from(key::ALLELE_COUNT))
    ///     .add_info(key::TOTAL_DEPTH, Map::<map::Info>::from(key::TOTAL_DEPTH))
    ///     .build();
    ///
    /// let string_maps = StringMaps::try_from(&header)?;
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    ///
    /// assert!(matches!(
    ///     info.get(&header, string_maps.strings(), key::ALLELE_COUNT),
    ///     Some(Ok(Some(Value::Integer(5))))
    /// ));
    ///
    /// assert!(info.get(&header, string_maps.strings(), &key::ANCESTRAL_ALLELE).is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn get<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_string_map: &'h StringStringMap,
        key: &str,
    ) -> Option<io::Result<Option<Value<'a>>>> {
        for result in self.iter(header, string_string_map) {
            match result {
                Ok((k, v)) => {
                    if k == key {
                        return Some(Ok(v));
                    }
                }
                Err(e) => return Some(Err(e)),
            }
        }

        None
    }

    /// Returns an iterator over all info fields.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_bcf::{header::StringMaps, record::Info};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map, Map},
    ///     record::info::field::key,
    ///     variant::record::info::field::Value,
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(key::ALLELE_COUNT, Map::<map::Info>::from(key::ALLELE_COUNT))
    ///     .add_info(key::TOTAL_DEPTH, Map::<map::Info>::from(key::TOTAL_DEPTH))
    ///     .build();
    ///
    /// let string_maps = StringMaps::try_from(&header)?;
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    /// let mut fields = info.iter(&header, string_maps.strings());
    ///
    /// assert!(matches!(
    ///     fields.next(),
    ///     Some(Ok((key::ALLELE_COUNT, Some(Value::Integer(5)))))
    /// ));
    ///
    /// assert!(matches!(
    ///     fields.next(),
    ///     Some(Ok((key::TOTAL_DEPTH, Some(Value::Integer(8)))))
    /// ));
    ///
    /// assert!(fields.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn iter<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_string_map: &'h StringStringMap,
    ) -> impl Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a {
        let mut reader = &self.buf[..];

        (0..self.len()).map(move |_| {
            read_field(&mut reader, header, string_string_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }

    /// Returns an iterator over all info values.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::{header::StringMaps, record::Info};
    /// use noodles_vcf::{
    ///     self as vcf,
    ///     header::record::value::{map, Map},
    ///     record::info::field::key,
    ///     variant::record::info::field::Value,
    /// };
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(key::ALLELE_COUNT, Map::<map::Info>::from(key::ALLELE_COUNT))
    ///     .add_info(key::TOTAL_DEPTH, Map::<map::Info>::from(key::TOTAL_DEPTH))
    ///     .build();
    ///
    /// let string_maps = StringMaps::try_from(&header)?;
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    ///
    /// let mut fields = info.values(&header, string_maps.strings());
    /// assert!(matches!(fields.next(), Some(Ok(Some(Value::Integer(5))))));
    /// assert!(matches!(fields.next(), Some(Ok(Some(Value::Integer(8))))));
    /// assert!(fields.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn values<'a, 'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_string_map: &'h StringStringMap,
    ) -> impl Iterator<Item = io::Result<Option<Value<'a>>>> + 'a {
        self.iter(header, string_string_map)
            .map(|result| result.map(|(_, value)| value))
    }

    pub(crate) fn set_field_count(&mut self, field_count: usize) {
        self.field_count = field_count;
    }
}

impl AsRef<[u8]> for Info {
    fn as_ref(&self) -> &[u8] {
        &self.buf
    }
}

impl AsMut<Vec<u8>> for Info {
    fn as_mut(&mut self) -> &mut Vec<u8> {
        &mut self.buf
    }
}
