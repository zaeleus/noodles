//! Indexed sequence reading with calculated newline skipping.

use std::io::{self, Read, Seek, SeekFrom};

use crate::fai;

/// Reads a sequence using raw byte calculation (no line parsing).
///
/// This reads all bytes in one syscall, then strips newlines in memory.
/// Much faster than line-by-line parsing for large sequences.
///
/// # Arguments
/// * `reader` - A seekable reader
/// * `index_record` - FAI record with offset/line metadata
/// * `start` - 0-based start position in the sequence
/// * `len` - Number of bases to read
pub fn read_sequence<R>(
    reader: &mut R,
    index_record: &fai::Record,
    start: u64,
    len: u64,
) -> io::Result<Vec<u8>>
where
    R: Read + Seek,
{
    if len == 0 {
        return Ok(Vec::new());
    }

    let line_bases = index_record.line_bases();
    let line_width = index_record.line_width();
    let terminator_len = (line_width - line_bases) as usize;

    // Calculate byte offset for start position
    let byte_offset =
        index_record.offset() + (start / line_bases) * line_width + (start % line_bases);

    // Calculate how many bytes to read (including newlines)
    let start_offset_in_line = (start % line_bases) as usize;
    let end = start + len;
    let end_offset_in_line = (end % line_bases) as usize;

    let start_line = start / line_bases;
    let end_line = if end_offset_in_line == 0 && end > start {
        (end / line_bases) - 1
    } else {
        end / line_bases
    };

    let total_bytes = if start_line == end_line {
        // All bases on same line
        len as usize
    } else {
        // First partial line + middle complete lines + last partial line
        let first_line_bases = (line_bases as usize) - start_offset_in_line;
        let middle_lines = (end_line - start_line - 1) as usize;
        let last_line_bases = if end_offset_in_line == 0 {
            line_bases as usize
        } else {
            end_offset_in_line
        };
        let num_terminators = (end_line - start_line) as usize;

        first_line_bases
            + (middle_lines * line_bases as usize)
            + last_line_bases
            + (num_terminators * terminator_len)
    };

    // One seek + one read
    reader.seek(SeekFrom::Start(byte_offset))?;
    let mut raw_bytes = vec![0u8; total_bytes];
    reader.read_exact(&mut raw_bytes)?;

    // Fast path: single line read has no newlines to strip
    if start_line == end_line {
        return Ok(raw_bytes);
    }

    // Strip newlines in memory
    let mut sequence = Vec::with_capacity(len as usize);
    let line_bases = line_bases as usize;

    let mut pos = 0;
    let mut offset_in_line = start_offset_in_line;

    while sequence.len() < len as usize && pos < raw_bytes.len() {
        // How many bases until end of current line?
        let bases_remaining_in_line = line_bases - offset_in_line;
        let bases_needed = (len as usize) - sequence.len();
        let bases_to_copy = bases_remaining_in_line
            .min(bases_needed)
            .min(raw_bytes.len() - pos);

        sequence.extend_from_slice(&raw_bytes[pos..pos + bases_to_copy]);
        pos += bases_to_copy;
        offset_in_line += bases_to_copy;

        // Skip terminator if we're at end of line and need more bases
        if offset_in_line == line_bases && sequence.len() < len as usize && pos < raw_bytes.len() {
            pos += terminator_len;
            offset_in_line = 0;
        }
    }

    Ok(sequence)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_record(
        name: &str,
        length: u64,
        offset: u64,
        line_bases: u64,
        line_width: u64,
    ) -> fai::Record {
        fai::Record::new(name, length, offset, line_bases, line_width)
    }

    #[test]
    fn test_read_full_sequence() -> io::Result<()> {
        // >sq0
        // ACGT
        // NNNN
        // AA
        let data = b">sq0\nACGT\nNNNN\nAA\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 10, 5, 4, 5);

        let seq = read_sequence(&mut reader, &record, 0, 10)?;
        assert_eq!(seq, b"ACGTNNNNAA");
        Ok(())
    }

    #[test]
    fn test_read_partial_middle() -> io::Result<()> {
        let data = b">sq0\nACGT\nNNNN\nAA\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 10, 5, 4, 5);

        // Read bases 2-5 (0-indexed): GT, NN
        let seq = read_sequence(&mut reader, &record, 2, 4)?;
        assert_eq!(seq, b"GTNN");
        Ok(())
    }

    #[test]
    fn test_read_spanning_lines() -> io::Result<()> {
        let data = b">sq0\nACGT\nNNNN\nAA\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 10, 5, 4, 5);

        // Read bases 3-6 (0-indexed): T, NNNN[0:3]
        let seq = read_sequence(&mut reader, &record, 3, 4)?;
        assert_eq!(seq, b"TNNN");
        Ok(())
    }

    #[test]
    fn test_read_single_base() -> io::Result<()> {
        let data = b">sq0\nACGT\nNNNN\nAA\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 10, 5, 4, 5);

        let seq = read_sequence(&mut reader, &record, 5, 1)?;
        assert_eq!(seq, b"N");
        Ok(())
    }

    #[test]
    fn test_read_crlf() -> io::Result<()> {
        // CRLF line endings
        let data = b">sq0\r\nACGT\r\nNNNN\r\nAA\r\n";
        let mut reader = Cursor::new(&data[..]);
        // offset 6 (">sq0\r\n" = 6 bytes, then sequence starts)
        let record = make_record("sq0", 10, 6, 4, 6);

        let seq = read_sequence(&mut reader, &record, 0, 10)?;
        assert_eq!(seq, b"ACGTNNNNAA");
        Ok(())
    }

    #[test]
    fn test_read_empty() -> io::Result<()> {
        let data = b">sq0\nACGT\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 4, 5, 4, 5);

        let seq = read_sequence(&mut reader, &record, 0, 0)?;
        assert_eq!(seq, b"");
        Ok(())
    }

    #[test]
    fn test_read_exact_line() -> io::Result<()> {
        let data = b">sq0\nACGT\nNNNN\n";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 8, 5, 4, 5);

        // Read exactly first line
        let seq = read_sequence(&mut reader, &record, 0, 4)?;
        assert_eq!(seq, b"ACGT");

        // Read exactly second line
        let mut reader = Cursor::new(&data[..]);
        let seq = read_sequence(&mut reader, &record, 4, 4)?;
        assert_eq!(seq, b"NNNN");
        Ok(())
    }

    #[test]
    fn test_read_last_partial_line() -> io::Result<()> {
        // Last line has only 2 bases (no trailing newline in some files)
        let data = b">sq0\nACGT\nNN";
        let mut reader = Cursor::new(&data[..]);
        let record = make_record("sq0", 6, 5, 4, 5);

        let seq = read_sequence(&mut reader, &record, 4, 2)?;
        assert_eq!(seq, b"NN");
        Ok(())
    }
}
