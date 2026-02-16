mod container;
mod file_id;
mod format_version;
mod magic_number;

use std::io::{self, Write};

use noodles_fasta as fasta;
use noodles_sam::{self as sam, header::ReferenceSequences};

use self::{
    container::write_container, file_id::write_file_id, format_version::write_format_version,
    magic_number::write_magic_number,
};
use crate::{FileDefinition, calculate_normalized_sequence_digest, file_definition::Version};

pub fn write_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    file_definition: &FileDefinition,
    header: &sam::Header,
    reference_required: bool,
) -> io::Result<()>
where
    W: Write,
{
    write_file_definition(writer, file_definition)?;
    write_file_header(
        writer,
        reference_sequence_repository,
        header,
        file_definition.version(),
        reference_required,
    )?;
    Ok(())
}

pub fn write_file_definition<W>(writer: &mut W, file_definition: &FileDefinition) -> io::Result<()>
where
    W: Write,
{
    write_magic_number(writer)?;
    write_format_version(writer, file_definition.version())?;
    write_file_id(writer, file_definition.file_id())?;
    Ok(())
}

pub fn write_file_header<W>(
    writer: &mut W,
    reference_sequence_repository: &fasta::Repository,
    header: &sam::Header,
    version: Version,
    reference_required: bool,
) -> io::Result<()>
where
    W: Write,
{
    let mut header = header.clone();

    add_missing_reference_sequence_checksums(
        reference_sequence_repository,
        header.reference_sequences_mut(),
        reference_required,
    )?;

    write_container(writer, &header, version)
}

pub(crate) fn add_missing_reference_sequence_checksums(
    reference_sequence_repository: &fasta::Repository,
    reference_sequences: &mut ReferenceSequences,
    reference_required: bool,
) -> io::Result<()> {
    use indexmap::map::Entry;
    use noodles_sam::header::record::value::map::reference_sequence::{Md5Checksum, tag};

    for (name, reference_sequence) in reference_sequences {
        if let Entry::Vacant(entry) = reference_sequence
            .other_fields_mut()
            .entry(tag::MD5_CHECKSUM)
        {
            match reference_sequence_repository.get(name).transpose()? {
                Some(sequence) => {
                    let checksum = calculate_normalized_sequence_digest(&sequence[..]);
                    entry.insert(Md5Checksum::from(checksum).to_string().into());
                }
                None if !reference_required => {
                    let checksum = [0u8; 16];
                    entry.insert(Md5Checksum::from(checksum).to_string().into());
                }
                None => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "missing reference sequence: {}",
                            String::from_utf8_lossy(name)
                        ),
                    ));
                }
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use bstr::BString;

    use super::*;

    #[test]
    fn test_add_missing_reference_sequence_checksums() -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use fasta::record::{Definition, Sequence};
        use sam::header::record::value::{
            Map,
            map::{ReferenceSequence, reference_sequence::tag},
        };

        let reference_sequences = vec![
            fasta::Record::new(
                Definition::new("sq0", None),
                Sequence::from(b"TTCACCCA".to_vec()),
            ),
            fasta::Record::new(
                Definition::new("sq1", None),
                Sequence::from(b"GATCTTACTTTTT".to_vec()),
            ),
        ];

        let sq0_md5_checksum = BString::from("be19336b7e15968f7ac7dc82493d9cd8");
        let sq1_md5_checksum = BString::from("d80f22a19aeeb623b3e4f746c762f21d");

        let repository = fasta::Repository::new(reference_sequences);

        let mut header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .add_reference_sequence(
                "sq1",
                Map::<ReferenceSequence>::builder()
                    .set_length(const { NonZero::new(13).unwrap() })
                    .insert(tag::MD5_CHECKSUM, sq1_md5_checksum.clone())
                    .build()?,
            )
            .build();

        add_missing_reference_sequence_checksums(
            &repository,
            header.reference_sequences_mut(),
            true,
        )?;

        let sq0 = header.reference_sequences().get(&b"sq0"[..]);
        assert_eq!(
            sq0.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
            Some(&sq0_md5_checksum)
        );

        let sq1 = header.reference_sequences().get(&b"sq1"[..]);
        assert_eq!(
            sq1.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
            Some(&sq1_md5_checksum)
        );

        Ok(())
    }

    #[test]
    fn test_add_missing_reference_sequence_checksums_without_reference()
    -> Result<(), Box<dyn std::error::Error>> {
        use std::num::NonZero;

        use sam::header::record::value::{
            Map,
            map::{ReferenceSequence, reference_sequence::tag},
        };

        let repository = fasta::Repository::default();

        let mut header = sam::Header::builder()
            .add_reference_sequence(
                "sq0",
                Map::<ReferenceSequence>::new(const { NonZero::new(8).unwrap() }),
            )
            .build();

        // reference_required=true should fail when reference is missing
        let result = add_missing_reference_sequence_checksums(
            &repository,
            header.reference_sequences_mut(),
            true,
        );
        assert!(result.is_err());

        // Ensure the MD5 field is vacant for the next call
        let sq0 = header
            .reference_sequences_mut()
            .get_mut(&b"sq0"[..])
            .unwrap();
        sq0.other_fields_mut().shift_remove(&tag::MD5_CHECKSUM);

        // reference_required=false should succeed and write a zero checksum
        add_missing_reference_sequence_checksums(
            &repository,
            header.reference_sequences_mut(),
            false,
        )?;

        let zero_md5 = BString::from("00000000000000000000000000000000");
        let sq0 = header.reference_sequences().get(&b"sq0"[..]);
        assert_eq!(
            sq0.and_then(|rs| rs.other_fields().get(&tag::MD5_CHECKSUM)),
            Some(&zero_md5)
        );

        Ok(())
    }
}
