#[derive(Clone, Debug)]
pub enum Feature {
    Bases(i32, Vec<u8>),
    Scores(i32, Vec<u8>),
    ReadBase(i32, u8, u8),
    Substitution(i32, u8),
    Insertion(i32, Vec<u8>),
    Deletion(i32, i32),
    InsertBase(i32, u8),
    QualityScore(i32, u8),
    ReferenceSkip(i32, i32),
    SoftClip(i32, Vec<u8>),
    Padding(i32, i32),
    HardClip(i32, i32),
}
