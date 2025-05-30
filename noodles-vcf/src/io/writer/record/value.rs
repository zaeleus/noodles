use percent_encoding::{AsciiSet, CONTROLS, PercentEncode, utf8_percent_encode};

// § 1.2 "Character encoding, non-printable characters and characters with special meaning" (2023-08-23)
const PERCENT_ENCODE_SET: &AsciiSet = &CONTROLS
    .add(b':')
    .add(b';')
    .add(b'=')
    .add(b'%')
    .add(b',')
    .add(b'\r')
    .add(b'\n')
    .add(b'\t');

pub(super) fn percent_encode(s: &str) -> PercentEncode<'_> {
    utf8_percent_encode(s, PERCENT_ENCODE_SET)
}
