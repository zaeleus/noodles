#[macro_use] extern crate clap;
extern crate noodles;

use clap::{App, Arg, SubCommand};
use noodles::convert;

fn main() {
    let convert_cmd = SubCommand::with_name("convert")
        .about("Converts between file formats")
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .value_name("file")
            .help("Output pathname")
            .required(true))
        .arg(Arg::with_name("input")
             .help("Input pathname")
             .required(true)
             .index(1));

    let matches = App::new(crate_name!())
        .version(crate_version!())
        .subcommand(convert_cmd)
        .get_matches();

    if let Some(matches) = matches.subcommand_matches("convert") {
        let src = matches.value_of("input").unwrap();
        let dst = matches.value_of("output").unwrap();
        convert(src, dst).unwrap();
    }
}
