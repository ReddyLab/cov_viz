use crate::data_structures::Options;
use std::collections::HashMap;

const DATABASE_URL_KEY: &str = "DATABASE_URL";

pub fn read_args<'a>(args: &'a Vec<String>, env_args: &'a HashMap<String, String>) -> Options<'a> {
    Options {
        output_location: &args[1],
        experiment_accession_id: &args[2],
        assembly_name: &args[3],
        connection_string: env_args.get(DATABASE_URL_KEY).unwrap(),
        bucket_size: match args.get(5) {
            Some(size) => u32::from_str_radix(size, 10).unwrap(),
            None => 2_000_000,
        },

        chromo: args.get(6),
    }
}
