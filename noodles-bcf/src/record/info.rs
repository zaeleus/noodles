use std::io;

use noodles_vcf as vcf;

use crate::header::StringMap;

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
    /// let string_map = bcf::header::StringMap::default();
    ///
    /// let vcf_info = bcf_info.try_into_vcf_record_info(&header, &string_map)?;
    /// assert!(vcf_info.is_empty());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn try_into_vcf_record_info(
        &self,
        header: &vcf::Header,
        string_map: &StringMap,
    ) -> io::Result<vcf::record::Info> {
        use crate::reader::record::site::read_info;
        let mut reader = &self.buf[..];
        read_info(&mut reader, header.infos(), string_map, self.len())
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

    /// Returns the field with the given key.
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_bcf::{header::StringMap, record::Info};
    /// use noodles_vcf::{self as vcf, record::info::{field::{Key, Value}, Field}, Header};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(vcf::header::Info::from(Key::AlleleCount))
    ///     .add_info(vcf::header::Info::from(Key::TotalDepth))
    ///     .build();
    ///
    /// let string_map = StringMap::from(&header);
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    ///
    /// assert_eq!(
    ///     info.get(&header, &string_map, &Key::AlleleCount).transpose()?,
    ///     Some(Field::new(Key::AlleleCount, Value::Integer(5)))
    /// );
    ///
    /// assert!(info.get(&header, &string_map, &Key::AncestralAllele).is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn get(
        &self,
        header: &vcf::Header,
        string_map: &StringMap,
        key: &vcf::record::info::field::Key,
    ) -> Option<io::Result<vcf::record::info::Field>> {
        for result in self.values(header, string_map) {
            match result {
                Ok(field) => {
                    if field.key() == key {
                        return Some(Ok(field));
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
    /// # use std::io;
    /// use noodles_bcf::{header::StringMap, record::Info};
    /// use noodles_vcf::{self as vcf, record::info::{field::{Key, Value}, Field}, Header};
    ///
    /// let header = vcf::Header::builder()
    ///     .add_info(vcf::header::Info::from(Key::AlleleCount))
    ///     .add_info(vcf::header::Info::from(Key::TotalDepth))
    ///     .build();
    ///
    /// let string_map = StringMap::from(&header);
    ///
    /// let data = vec![
    ///     0x11, 0x01, 0x11, 0x05, // AC=5
    ///     0x11, 0x02, 0x11, 0x08, // DP=8
    /// ];
    ///
    /// let info = Info::new(data, 2);
    ///
    /// let mut fields = info.values(&header, &string_map);
    ///
    /// assert_eq!(fields.next().transpose()?, Some(Field::new(Key::AlleleCount, Value::Integer(5))));
    /// assert_eq!(fields.next().transpose()?, Some(Field::new(Key::TotalDepth, Value::Integer(8))));
    /// assert!(fields.next().is_none());
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn values<'a>(
        &'a self,
        header: &'a vcf::Header,
        string_map: &'a StringMap,
    ) -> impl Iterator<Item = io::Result<vcf::record::info::Field>> + 'a {
        use crate::reader::record::site::info::read_info_field;
        let mut reader = &self.buf[..];
        (0..self.len()).map(move |_| read_info_field(&mut reader, header.infos(), string_map))
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
