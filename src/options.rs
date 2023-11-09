use std::collections::HashMap;
use std::env;
use std::path::PathBuf;

const DATABASE_URL_KEY: &str = "DATABASE_URL";

#[derive(Debug)]
pub struct Options {
    pub cov_output_location: PathBuf,
    pub features_output_location: PathBuf,
    pub analysis_accession_id: String,
    pub assembly_name: String,
    pub connection_string: String,
    pub bucket_size: u32,
    pub chromo: Option<String>,
}

impl Options {
    pub fn get() -> Self {
        let env_args: HashMap<String, String> = env::vars().collect();
        let args: Vec<String> = env::args().collect();

        let output_location = &args[1];
        let chromo = match args.get(5) {
            Some(chrom_name) => Some(chrom_name.to_string()),
            None => None,
        };

        let cov_output_file = match &chromo {
            Some(chrom_name) => format!("level2_{}.ecd", chrom_name.strip_prefix("chr").unwrap()),
            None => "level1.ecd".to_string(),
        };
        let cov_path: PathBuf = [output_location, &cov_output_file].iter().collect();

        let feat_output_file = match &chromo {
            Some(chrom_name) => format!("level2_{}.fd", chrom_name.strip_prefix("chr").unwrap()),
            None => "level1.fd".to_string(),
        };
        let feat_path: PathBuf = [output_location, &feat_output_file].iter().collect();

        Options {
            cov_output_location: cov_path,
            features_output_location: feat_path,
            analysis_accession_id: args[2].clone(),
            assembly_name: args[3].clone(),
            bucket_size: match args.get(4) {
                Some(size) => u32::from_str_radix(size, 10).unwrap(),
                None => 2_000_000,
            },
            chromo: chromo,
            connection_string: env_args.get(DATABASE_URL_KEY).unwrap().to_string(),
        }
    }
}
