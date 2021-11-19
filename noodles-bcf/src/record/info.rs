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
