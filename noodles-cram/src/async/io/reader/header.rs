pub(super) mod container;
mod file_id;
mod format_version;
mod magic_number;

use tokio::io::{self, AsyncRead};

use self::{
    file_id::read_file_id, format_version::read_format_version, magic_number::read_magic_number,
};
use crate::FileDefinition;

pub(super) async fn read_file_definition<R>(reader: &mut R) -> io::Result<FileDefinition>
where
    R: AsyncRead + Unpin,
{
    use crate::io::reader::header::magic_number;

    read_magic_number(reader)
        .await
        .and_then(magic_number::validate)?;

    let version = read_format_version(reader).await?;
    let file_id = read_file_id(reader).await?;

    Ok(FileDefinition::new(version, file_id))
}
