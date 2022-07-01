use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use bincode::Options as BincodeOptions;

use crate::data_structures::CoverageData;

pub fn serialize(data: &CoverageData, chromo: Option<&String>, output_location: &String) {
    let bincode_options = bincode::DefaultOptions::new().with_no_limit();
    let output_file = match chromo {
        Some(chrom_name) => format!("level2_{}.bin", chrom_name.strip_prefix("chr").unwrap()),
        None => "level1.bin".to_string(),
    };
    let path: PathBuf = [output_location, &output_file].iter().collect();

    let mut writer = BufWriter::new(File::create(path).unwrap());
    let encoded_start_time = Instant::now();
    bincode_options.serialize_into(&mut writer, &data).unwrap();
    match writer.flush() {
        Ok(_) => (),
        Err(e) => println!("Error flushing during serialization: {:?}", e),
    };
    println!(
        "Serialization Finished... {:}ms",
        encoded_start_time.elapsed().as_millis()
    );
}
