mod build_data;
mod data_structures;
mod options;

use postgres::{Client, NoTls};

use crate::build_data::build_data;
use crate::options::Options;

type DbID = i64;

fn main() {
    let options = Options::get();

    let mut client = match Client::connect(&options.connection_string, NoTls) {
        Ok(client) => client,
        Err(e) => {
            panic!("{}", e)
        }
    };

    match build_data(&options, &mut client) {
        Ok(data) => data.serialize(&options.output_location),
        Err(e) => eprintln!("{}", e),
    };
}
