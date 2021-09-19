mod tag;
mod value;

use self::{tag::read_tag, value::read_value};

use tokio::io::{self, AsyncBufRead};

use crate::record::data::Field;

pub async fn read_field<R>(reader: &mut R) -> io::Result<Field>
where
    R: AsyncBufRead + Unpin,
{
    let tag = read_tag(reader).await?;
    let ty = value::read_type(reader).await?;
    let value = read_value(reader, ty).await?;
    Ok(Field::new(tag, value))
}
