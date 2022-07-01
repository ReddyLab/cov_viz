use serde::{Deserialize, Serialize};

use crate::data_structures::Interval;

#[derive(Serialize, Deserialize, Debug)]
pub struct ChromosomeData {
    pub chrom: String,
    pub bucket_size: u32,
    pub target_intervals: Vec<Interval>,
    pub source_intervals: Vec<Interval>,
}

impl ChromosomeData {
    pub fn from(chrom: &str, bucket_size: u32) -> Self {
        let cd = ChromosomeData {
            chrom: chrom.to_string(),
            bucket_size: bucket_size,
            source_intervals: Vec::new(),
            target_intervals: Vec::new(),
        };
        cd
    }

    // Uncomment if needed. This can be useful for debugging
    // pub fn print_info(&self) {
    //     println!(
    //         "\n{} {} {} {}",
    //         self.chrom,
    //         self.bucket_size,
    //         self.target_intervals.len(),
    //         self.source_intervals.len()
    //     );
    //     print!("    ");
    //     for interval in &self.source_intervals {
    //         let count = interval.values.borrow().facets.len();
    //         if count > 0 {
    //             print!("{} {} ", interval.start, count);
    //         }
    //     }
    //     print!("\n    ");
    //     for interval in &self.target_intervals {
    //         let count = interval.values.borrow().facets.len();
    //         if count > 0 {
    //             print!("{} {} ", interval.start, count);
    //         }
    //     }
    //     println!("");
    // }
}
