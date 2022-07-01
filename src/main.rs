mod args;
mod build_data;
mod data_structures;
mod serialize;
mod utils;

use std::env;
use std::collections::HashMap;

use postgres::{Client, Error, NoTls};

use crate::serialize::serialize;
use crate::build_data::build_data;
use crate::args::read_args;

type DbID = i64;

fn main() -> Result<(), Error> {
    let env_args: HashMap<String, String> = env::vars().collect();
    let args: Vec<String> = env::args().collect();
    let options = read_args(&args, &env_args);

    let mut client = Client::connect(options.connection_string,
        NoTls,
    )?;

    let data = build_data(&options, &mut client)?;

    serialize(&data, options.chromo, options.output_location);

    Ok(())
}
