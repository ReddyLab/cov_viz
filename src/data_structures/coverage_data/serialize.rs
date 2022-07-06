use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use bincode::Options as BincodeOptions;

use crate::data_structures::CoverageData;

impl CoverageData {
    pub fn serialize(&self, output_path: &PathBuf) {
        let bincode_options = bincode::DefaultOptions::new().with_no_limit();

        let mut writer = BufWriter::new(File::create(output_path).unwrap());
        bincode_options.serialize_into(&mut writer, self).unwrap();
        match writer.flush() {
            Ok(_) => (),
            Err(e) => eprintln!("Error flushing during serialization: {:?}", e),
        };
    }

    pub fn deserialize(file_path: &PathBuf) -> Result<Self, bincode::Error> {
        let bincode_options = bincode::DefaultOptions::new().with_no_limit();
        let raw_bytes = fs::read(file_path).unwrap();
        bincode_options.deserialize(&raw_bytes)
    }
}
