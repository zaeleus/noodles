//! FASTX record

/* std use */

/* crates use */

/* mod declaration */

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

    pub fn is_fastq(&self) -> bool {
        self.quality.is_some()
    }

    pub fn fastq2fasta(&mut self) {
        self.quality = None
    }
}
