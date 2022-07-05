mod args;
mod build_data;
mod data_structures;

use std::collections::HashMap;
use std::env;
use std::path::PathBuf;

use postgres::{Client, Error, NoTls};

use crate::args::read_args;
use crate::build_data::build_data;

type DbID = i64;

fn main() -> Result<(), Error> {
    let env_args: HashMap<String, String> = env::vars().collect();
    let args: Vec<String> = env::args().collect();
    let options = read_args(&args, &env_args);

    let mut client = Client::connect(options.connection_string, NoTls)?;

    let data = build_data(&options, &mut client)?;

    let output_file = match options.chromo {
        Some(chrom_name) => format!("level2_{}.bin", chrom_name.strip_prefix("chr").unwrap()),
        None => "level1.bin".to_string(),
    };
    let path: PathBuf = [options.output_location, &output_file].iter().collect();

    data.serialize(&path);

    Ok(())
}
