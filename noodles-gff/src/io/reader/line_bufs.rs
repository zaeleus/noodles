use std::io::{self, BufRead};

use super::Reader;
use crate::{Line, LineBuf, feature::RecordBuf, line::Kind};

/// An iterator over lines of a GFF reader.
///
/// When using this, the caller is responsible to stop reading at either EOF or when the `FASTA`
/// directive is read, whichever comes first.
///
/// This is created by calling [`Reader::lines`].
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
                Kind::Comment => Some(Ok(LineBuf::Comment(self.line.as_ref().into()))),
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
