mod block;
mod header;

pub use self::{block::write_block, header::write_header};

use std::io::{self, Write};

use crate::{container::Block, Container};

pub fn write_container<W>(writer: &mut W, container: &Container) -> io::Result<()>
where
    W: Write,
{
    write_header(writer, container.header())?;
    write_blocks(writer, container.blocks())?;
    Ok(())
}

fn write_blocks<W>(writer: &mut W, blocks: &[Block]) -> io::Result<()>
where
    W: Write,
{
    for block in blocks {
        write_block(writer, block)?;
    }

    Ok(())
}
