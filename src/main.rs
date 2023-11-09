mod build_data;
mod options;

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
        Ok((coverage, features)) => {
            coverage.serialize(&options.cov_output_location);
            features.serialize(&options.features_output_location);
        }
        Err(e) => eprintln!("{}", e),
    };
}
