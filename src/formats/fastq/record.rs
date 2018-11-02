#[derive(Default, Debug)]
pub struct Record {
    name: Vec<u8>,
    sequence: Vec<u8>,
    plus_line: Vec<u8>,
    quality: Vec<u8>,
}

impl Record {
    pub fn new() -> Record {
        Record::default()
    }

    pub fn name(&self) -> &[u8] {
        &self.name
    }

    pub fn name_mut(&mut self) -> &mut Vec<u8> {
        &mut self.name
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn sequence_mut(&mut self) -> &mut Vec<u8> {
        &mut self.sequence
    }

    pub fn plus_line(&self) -> &[u8] {
        &self.plus_line
    }

    pub fn plus_line_mut(&mut self) -> &mut Vec<u8> {
        &mut self.plus_line
    }

    pub fn quality(&self) -> &[u8] {
        &self.quality
    }

    pub fn quality_mut(&mut self) -> &mut Vec<u8> {
        &mut self.quality
    }

    pub fn clear(&mut self) {
        self.name.clear();
        self.sequence.clear();
        self.plus_line.clear();
        self.quality.clear();
    }
}
