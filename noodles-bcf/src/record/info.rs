mod field;

use std::io;

use noodles_vcf::{self as vcf, variant::record::info::field::Value};

use self::field::read_field;
use crate::header::string_maps::StringStringMap;

/// BCF record info.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Info<'a> {
    src: &'a [u8],
    field_count: usize,
}

impl<'a> Info<'a> {
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
    pub fn try_into_vcf_record_info<'h: 'a>(
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

    pub(super) fn new(src: &'a [u8], field_count: usize) -> Self {
        Self { src, field_count }
    }

    /// Returns whether there are any info fields.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns the number of info fields.
    pub fn len(&self) -> usize {
        self.field_count
    }

    /// Returns the value with the given key.
    pub fn get<'h: 'a>(
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

    /// Returns an iterator over info fields.
    pub fn iter<'h: 'a>(
        &'a self,
        header: &'h vcf::Header,
        string_string_map: &'h StringStringMap,
    ) -> impl Iterator<Item = io::Result<(&'a str, Option<Value<'a>>)>> + 'a {
        let mut src = self.as_ref();

        (0..self.len()).map(move |_| {
            read_field(&mut src, header, string_string_map)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
        })
    }
}

impl<'a> AsRef<[u8]> for Info<'a> {
    fn as_ref(&self) -> &[u8] {
        self.src
    }
}
