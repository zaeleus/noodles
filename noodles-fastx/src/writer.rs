//! Define a fastx writer

/* std use */

/* crates use */

/* project use */
use crate::record::Record;

/// A FASTX writer
pub struct Writer<W> {
    inner: W,
}

impl<W> Writer<W>
where
    W: std::io::Write,
{
    /// Creates a FASTX writer.
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    /// Returns a reference to the underlying writer.
    ///
    /// # Examples
    ///
    /// ```
    /// use noodles_fastx as fastx;
    /// let writer = fastx::Writer::new(Vec::new());
    /// assert!(writer.get_ref().is_empty());
    /// ```
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    /// Writes a FASTX record.
    pub fn write_record(&mut self, record: &Record) -> std::io::Result<usize> {
        record.to_writer(&mut self.inner)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_record() -> std::io::Result<()> {
        let mut writer = Writer::new(Vec::new());

        let mut record = Record::default();
        record.name_mut().extend(b"1");
        record.sequence_mut().extend(b"ACTCA");
        writer.write_record(&record)?;

        let mut record = Record::default();
        record.name_mut().extend(b"2");
        record.sequence_mut().extend(b"ACTCA");
        *record.description_mut() = Some(b"second record".to_vec());
        *record.second_description_mut() = Some(b"second description".to_vec());
        *record.quality_mut() = Some(b"!!;!!".to_vec());
        writer.write_record(&record)?;

        let mut record = Record::default();
        record.name_mut().extend(b"3");
        record.sequence_mut().extend(b"ACTCA");
        *record.description_mut() = Some(b"second record".to_vec());
        *record.second_description_mut() = Some(b"second description".to_vec());
        writer.write_record(&record)?;

        let mut record = Record::default();
        record.name_mut().extend(b"4");
        record.sequence_mut().extend(b"ACTCA");
        *record.description_mut() = Some(b"second record".to_vec());
        *record.quality_mut() = Some(b"!!;!!".to_vec());
        writer.write_record(&record)?;

        let mut record = Record::default();
        record.name_mut().extend(b"5");
        record.sequence_mut().extend(b"ACTCA");
        *record.second_description_mut() = Some(b"second description".to_vec());
        *record.quality_mut() = Some(b"!!;!!".to_vec());
        writer.write_record(&record)?;

        let expected = b"\
>1
ACTCA
@2 second record
ACTCA
+second description
!!;!!
>3 second record
ACTCA
@4 second record
ACTCA
+
!!;!!
@5
ACTCA
+second description
!!;!!
"
        .to_vec();

        assert_eq!(
            String::from_utf8(expected),
            String::from_utf8(writer.get_ref().to_vec())
        );

        Ok(())
    }
}
