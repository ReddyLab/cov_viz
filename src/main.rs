mod build_data;
mod data_structures;
mod options;

use std::time::Instant;

use data_structures::CoverageData;
use postgres::{Client, NoTls};

use crate::build_data::build_data;
use crate::options::Options;

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

    let now = Instant::now();
    let _ = CoverageData::deserialize(&options.output_location);
    println!("Time to decode: {}ms", now.elapsed().as_millis());
}
