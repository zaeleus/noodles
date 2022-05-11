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

pub fn write_eof_container<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    // § 9 "End of file container" (2022-04-12)
    static EOF: [u8; 38] = [
        0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00,
        0x00, 0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00,
        0x01, 0x00, 0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
    ];

    writer.write_all(&EOF)
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
