use std::io::{self, BufRead};

use super::Reader;
use crate::{Line, LineBuf, feature::RecordBuf, line::Kind};

/// An iterator over lines of a GFF reader.
///
/// When using this, the caller is responsible to stop reading at either EOF or when the `FASTA`
/// directive is read, whichever comes first.
///
/// This is created by calling [`Reader::line_bufs`].
pub struct LineBufs<'a, R> {
    inner: &'a mut Reader<R>,
    line: Line,
}

impl<'a, R> LineBufs<'a, R>
where
    R: BufRead,
{
    pub(crate) fn new(inner: &'a mut Reader<R>) -> Self {
        Self {
            inner,
            line: Line::default(),
        }
    }
}

impl<R> Iterator for LineBufs<'_, R>
where
    R: BufRead,
{
    type Item = io::Result<LineBuf>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.read_line(&mut self.line) {
            Ok(0) => None,
            Ok(_) => match self.line.kind() {
                Kind::Directive => {
                    // SAFETY: `self.line` is a directive.
                    let directive = self.line.as_directive().unwrap();
                    Some(Ok(LineBuf::Directive(directive.into())))
                }
                Kind::Comment => {
                    // SAFETY: `self.line` is a comment.
                    let comment = self.line.as_comment().unwrap();
                    Some(Ok(LineBuf::Comment(comment.into())))
                }
                Kind::Record => Some(
                    self.line
                        .as_record()
                        .unwrap() // SAFETY: `self.line` is a record.
                        .and_then(|record| {
                            RecordBuf::try_from_feature_record(&record).map(LineBuf::Record)
                        }),
                ),
            },
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;
    use crate::{DirectiveBuf, directive_buf::Value};

    #[test]
    fn test_next() -> io::Result<()> {
        let mut reader =
            Reader::new(&b"##gff-version 3\n#noodles\n.\t.\t.\t1\t1\t.\t.\t.\t.\n"[..]);

        let iter = LineBufs::new(&mut reader);
        let actual: Vec<_> = iter.collect::<io::Result<_>>()?;

        let expected = [
            LineBuf::Directive(DirectiveBuf::new(
                "gff-version",
                Some(Value::String(BString::from("3"))),
            )),
            LineBuf::Comment(BString::from("noodles")),
            LineBuf::Record(RecordBuf::default()),
        ];

        assert_eq!(actual, expected);

        Ok(())
    }
}
