mod build_data;
mod data_structures;
mod options;

use postgres::{Client, Error, NoTls};

use crate::build_data::build_data;
use crate::options::Options;

type DbID = i64;

fn main() -> Result<(), Error> {
    let options = Options::get();

    let mut client = Client::connect(&options.connection_string, NoTls)?;

    let data = build_data(&options, &mut client)?;

    data.serialize(&options.output_location);

    Ok(())
}
