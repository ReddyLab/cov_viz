pub struct Options<'a> {
    pub output_location: &'a String,
    pub experiment_accession_id: &'a String,
    pub assembly_name: &'a String,
    pub connection_string: &'a String,
    pub bucket_size: u32,
    pub chromo: Option<&'a String>,
}
