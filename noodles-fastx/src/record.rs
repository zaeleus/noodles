//! FASTX record

/* std use */

/* crates use */
use bstr::ByteSlice;

/* mod declaration */

/// Error could be generate by FASTX reading and writing
pub enum Error {
    /// Not the good sequence prefix
    MissingPrefix,
    /// No id
    MissingName,
    /// No sequence
    MissingSequence,
    /// No second description (fastq only)
    MissingSecondDescription,
    /// No or empty quality string (fastq only)
    MissingQuality,
    /// Too short buffer
    TooShortBuffer,
}

// Convert an std::io::Error in Error is impossible this clippy warning isn't intrusting here
#[allow(clippy::from_over_into)]
impl Into<std::io::Error> for Error {
    fn into(self) -> std::io::Error {
        match self {
            Error::MissingPrefix => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The prefix ('>' or '@') is missing",
            ),
            Error::MissingSequence => {
                std::io::Error::new(std::io::ErrorKind::InvalidData, "The sequence is missing")
            }
            Error::MissingName => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The sequence name is missing",
            ),
            Error::MissingSecondDescription => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The second description is missing",
            ),
            Error::MissingQuality => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The quality string is missing",
            ),
            Error::TooShortBuffer => std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "The reading buffer is shorter than line",
            ),
        }
    }
}

/// Fastx record struct
#[derive(Clone, Debug, Eq, PartialEq, Default)]
pub struct Record {
    name: Vec<u8>,
    description: Option<Vec<u8>>,
    second_description: Option<Vec<u8>>,
    sequence: Vec<u8>,
    quality: Option<Vec<u8>>,
}

impl Record {
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    pub fn description(&self) -> Option<&[u8]> {
        self.description.as_deref()
    }

    pub fn second_description(&self) -> Option<&[u8]> {
        self.second_description.as_deref()
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn quality(&self) -> Option<&[u8]> {
        self.quality.as_deref()
    }

    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.name
    }

    pub fn description_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.description
    }

    pub fn second_description_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.second_description
    }

    pub fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    pub fn quality_mut(&mut self) -> &mut Option<Vec<u8>> {
        &mut self.quality
    }

    pub fn from_reader<R>(&mut self, reader: &mut R) -> std::io::Result<usize>
    where
        R: std::io::BufRead,
    {
        self.clear();

        let buffer = reader.fill_buf()?;

        if buffer.is_empty() {
            return Ok(0);
        }

        match buffer[0] {
            b'>' => read_fasta(self, reader),
            b'@' => read_fastq(self, reader),
            _ => Err(Error::MissingPrefix.into()),
        }
    }

    pub fn to_writer<W>(&self, writer: &mut W) -> std::io::Result<usize>
    where
        W: std::io::Write,
    {
        let mut buffer = Vec::new();

        match (
            self.description.as_ref(),
            self.second_description.as_ref(),
            self.quality.as_ref(),
        ) {
            (None, None, None) => {
                buffer.push(b'>');
                buffer.extend(&self.name);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
            }
            (None, None, Some(qual)) => {
                buffer.push(b'@');
                buffer.extend(&self.name);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
                buffer.extend(b"+\n");
                buffer.extend(qual);
                buffer.push(b'\n');
            }
            (None, Some(_), None) => {
                buffer.push(b'>');
                buffer.extend(&self.name);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
            }
            (None, Some(second_desc), Some(qual)) => {
                buffer.push(b'@');
                buffer.extend(&self.name);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
                buffer.push(b'+');
                buffer.extend(second_desc);
                buffer.push(b'\n');
                buffer.extend(qual);
                buffer.push(b'\n')
            }
            (Some(desc), None, None) => {
                buffer.push(b'>');
                buffer.extend(&self.name);
                buffer.push(b' ');
                buffer.extend(desc);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
            }
            (Some(desc), None, Some(qual)) => {
                buffer.push(b'@');
                buffer.extend(&self.name);
                buffer.push(b' ');
                buffer.extend(desc);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.extend(b"\n+\n");
                buffer.extend(qual);
                buffer.push(b'\n');
            }
            (Some(desc), Some(_), None) => {
                buffer.push(b'>');
                buffer.extend(&self.name);
                buffer.push(b' ');
                buffer.extend(desc);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
            }
            (Some(desc), Some(second_desc), Some(qual)) => {
                buffer.push(b'@');
                buffer.extend(&self.name);
                buffer.push(b' ');
                buffer.extend(desc);
                buffer.push(b'\n');
                buffer.extend(&self.sequence);
                buffer.push(b'\n');
                buffer.push(b'+');
                buffer.extend(second_desc);
                buffer.push(b'\n');
                buffer.extend(qual);
                buffer.push(b'\n');
            }
        }

        writer.write_all(&buffer[..])?;

        Ok(buffer.len())
    }

    pub fn is_fastq(&self) -> bool {
        self.quality.is_some()
    }

    pub fn fastq2fasta(&mut self) {
        self.quality = None
    }

    pub fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();

        self.description = None;
        self.second_description = None;
        self.quality = None;
    }
}

fn get_end_line(buffer: &[u8]) -> std::io::Result<(usize, usize)> {
    match buffer.find_byte(b'\n') {
        Some(pos) if buffer[..pos].ends_with(&[b'\r']) => Ok((pos - 1, pos + 1)),
        Some(pos) => Ok((pos, pos + 1)),
        None => Err(Error::TooShortBuffer.into()),
    }
}

fn read_fasta<R>(record: &mut Record, reader: &mut R) -> std::io::Result<usize>
where
    R: std::io::BufRead,
{
    let mut bytes_read = 0;

    let buffer = reader.fill_buf()?;

    // Found name and description
    let (eol, bon) = get_end_line(buffer)?;
    let line = &buffer[..eol];
    if buffer[..eol].len() <= 1 {
        return Err(Error::MissingName.into());
    }
    match line[..].find(&[b' ']) {
        Some(i) => {
            record.name_mut().extend(&line[1..i]);
            *record.description_mut() = Some(line[i + 1..].to_vec());
        }
        None => record.name_mut().extend(&line[1..]),
    }
    bytes_read += bon;
    reader.consume(bon);

    // Found sequence
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() || buffer[0] == b'>' || buffer[0] == b'@' {
            break;
        }

        let (eol, bon) = get_end_line(buffer)?;
        record.sequence_mut().extend(&buffer[..eol]);

        reader.consume(bon);
        bytes_read += bon
    }

    if record.sequence().is_empty() {
        return Err(Error::MissingSequence.into());
    }

    Ok(bytes_read)
}

fn read_fastq<R>(record: &mut Record, reader: &mut R) -> std::io::Result<usize>
where
    R: std::io::BufRead,
{
    let mut bytes_read = 0;

    let buffer = reader.fill_buf()?;

    // Found name and description
    let (eol, bon) = get_end_line(buffer)?;
    let line = &buffer[..eol];
    if buffer[..eol].len() <= 1 {
        return Err(Error::MissingName.into());
    }
    match line[..].find(&[b' ']) {
        Some(i) => {
            record.name_mut().extend(&line[1..i]);
            *record.description_mut() = Some(line[i + 1..].to_vec());
        }
        None => record.name_mut().extend(&line[1..]),
    }
    bytes_read += bon;
    reader.consume(bon);

    // Found sequence
    let buffer = reader.fill_buf()?;
    let (eol, bon) = get_end_line(buffer)?;
    let line = &buffer[..eol];
    if buffer[..eol].len() <= 2 {
        return Err(Error::MissingSequence.into());
    }
    record.sequence_mut().extend(&line[..eol]);
    bytes_read += bon;
    reader.consume(bon);

    // Found second description
    let buffer = reader.fill_buf()?;
    let (eol, bon) = get_end_line(buffer)?;
    let line = &buffer[..eol];
    if line[0] != b'+' {
        return Err(Error::MissingSecondDescription.into());
    }
    if line[1..].is_empty() {
        *record.second_description_mut() = None;
    } else {
        *record.second_description_mut() = Some(buffer[1..eol].to_vec());
    }
    bytes_read += bon;
    reader.consume(bon);

    // Found quality
    let buffer = reader.fill_buf()?;
    let (eol, bon) = get_end_line(buffer)?;
    let line = &buffer[..eol];
    if buffer[..eol].is_empty() {
        return Err(Error::MissingQuality.into());
    }
    *record.quality_mut() = Some(line[..eol].to_vec());
    bytes_read += bon;
    reader.consume(bon);

    Ok(bytes_read)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_record() -> std::io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATC
ATGA
>1
ACGT
";

        let mut reader = &data[..];
        let mut record = Record::default();

        assert_eq!(25, record.from_reader(&mut reader)?);
        assert_eq!(b"noodles:1/1", record.name());
        assert_eq!(None, record.description());
        assert_eq!(b"AGCT", record.sequence());
        assert_eq!(None, record.second_description());
        assert_eq!(Some(b"abcd".as_ref()), record.quality());

        assert_eq!(37, record.from_reader(&mut reader)?);
        assert_eq!(b"noodles:", record.name());
        assert_eq!(b"2/1", record.description().unwrap());
        assert_eq!(b"TCGA", record.sequence());
        assert_eq!(Some(b"dcba".as_ref()), record.quality());
        assert_eq!(b"noodles:2/1", record.second_description().unwrap());

        assert_eq!(48, record.from_reader(&mut reader)?);
        assert_eq!(b"fasta", record.name());
        assert_eq!(
            b"sequence in middle of my fastq",
            record.description().unwrap()
        );
        assert_eq!(b"GATCATGA", record.sequence());
        assert_eq!(None, record.quality());
        assert_eq!(None, record.second_description());

        assert_eq!(8, record.from_reader(&mut reader)?);
        assert_eq!(b"1", record.name());
        assert_eq!(None, record.description());
        assert_eq!(b"ACGT", record.sequence());
        assert_eq!(None, record.quality());
        assert_eq!(None, record.second_description());

        assert_eq!(0, record.from_reader(&mut reader)?);

        Ok(())
    }

    #[test]
    fn bad_record() -> std::io::Result<()> {
        let bad_record = b"\
!noodles
ACGT
+
abcd
";

        let mut reader = &bad_record[..];
        let mut record = Record::default();

        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The prefix ('>' or '@') is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@
ACGT
+
abcd
";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence name is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles

+
abcd
";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles
ACTG
abcd
";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The second description is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
@noodles
ACTG
+

";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The quality string is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
>missing sequence

";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
>
ACCTA
";

        let mut reader = &bad_record[..];
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The sequence name is missing"
                ))
            ),
            format!("{:?}", result.err())
        );

        let bad_record = b"\
>aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
";
        let mut reader = std::io::BufReader::with_capacity(10, &bad_record[..]);
        let result = record.from_reader(&mut reader);
        assert_eq!(
            format!(
                "{:?}",
                Some(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "The reading buffer is shorter than line"
                ))
            ),
            format!("{:?}", result.err())
        );

        Ok(())
    }

    #[test]
    fn read_write() -> std::io::Result<()> {
        let data = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATC
ATGA
>aue
AACAC
@2
AGCT
+strange description
!!!!!!!!!
";

        let mut reader = &data[..];
        let mut writer: Vec<u8> = Vec::new();
        let mut record = Record::default();

        record.from_reader(&mut reader)?;
        assert!(record.is_fastq());
        record.to_writer(&mut writer)?;

        record.from_reader(&mut reader)?;
        assert!(record.is_fastq());
        record.to_writer(&mut writer)?;

        record.from_reader(&mut reader)?;
        assert!(!record.is_fastq());
        record.to_writer(&mut writer)?;

        record.from_reader(&mut reader)?;
        assert!(!record.is_fastq());
        record.to_writer(&mut writer)?;

        record.from_reader(&mut reader)?;
        assert!(record.is_fastq());
        record.fastq2fasta();
        assert!(!record.is_fastq());
        record.to_writer(&mut writer)?;

        assert_eq!(0, record.from_reader(&mut reader).unwrap());

        let expected_result = b"\
@noodles:1/1
AGCT
+
abcd
@noodles: 2/1
TCGA
+noodles:2/1
dcba
>fasta sequence in middle of my fastq
GATCATGA
>aue
AACAC
>2
AGCT
";

        assert_eq!(
            String::from_utf8(expected_result.to_vec()).unwrap(),
            String::from_utf8(writer).unwrap()
        );

        Ok(())
    }
}
