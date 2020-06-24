use std::io::{self, BufRead};

const HEADER_PREFIX: u8 = b'#';
const NEWLINE: u8 = b'\n';

pub struct Reader<R> {
    inner: R,
}

impl<R> Reader<R>
where
    R: BufRead,
{
    pub fn new(inner: R) -> Self {
        Self { inner }
    }

    /// Reads the raw GFF header.
    ///
    /// The GFF header is all the comments (lines starting with `#`) preceding the list of records.
    ///
    /// The position of the stream is expected to be at the start.
    ///
    /// This returns the raw header as a [`String`].
    ///
    /// [`String`]: https://doc.rust-lang.org/std/string/struct.String.html
    ///
    /// # Examples
    ///
    /// ```
    /// # use std::io;
    /// use noodles_gff as gff;
    ///
    /// let data = b"##format: gtf
    /// sq0\tNOODLES\tgene\t8\t13\t.\t+\t.\tgene_id \"ndls0\"; gene_name \"g0\"
    /// ";
    /// let mut reader = gff::Reader::new(&data[..]);
    ///
    /// let header = reader.read_header()?;
    /// assert_eq!(header, "##format: gtf\n");
    /// # Ok::<_, io::Error>(())
    /// ```
    pub fn read_header(&mut self) -> io::Result<String> {
        let mut header_buf = Vec::new();
        let mut eol = false;

        loop {
            let buf = self.inner.fill_buf()?;

            if eol && !buf.is_empty() && buf[0] != HEADER_PREFIX {
                break;
            }

            let (read_eol, len) = match buf.iter().position(|&b| b == NEWLINE) {
                Some(i) => {
                    header_buf.extend(&buf[..=i]);
                    (true, i + 1)
                }
                None => {
                    header_buf.extend(buf);
                    (false, buf.len())
                }
            };

            eol = read_eol;
            self.inner.consume(len);
        }

        String::from_utf8(header_buf).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
    }
}
