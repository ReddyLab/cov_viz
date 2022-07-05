use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::time::Instant;

use bincode::Options as BincodeOptions;

use crate::data_structures::CoverageData;

impl CoverageData {
    pub fn serialize(&self, output_path: &PathBuf) {
        let bincode_options = bincode::DefaultOptions::new().with_no_limit();

        let mut writer = BufWriter::new(File::create(output_path).unwrap());
        let encoded_start_time = Instant::now();
        bincode_options.serialize_into(&mut writer, self).unwrap();
        match writer.flush() {
            Ok(_) => (),
            Err(e) => println!("Error flushing during serialization: {:?}", e),
        };
        println!(
            "Serialization Finished... {:}ms",
            encoded_start_time.elapsed().as_millis()
        );
    }
}
