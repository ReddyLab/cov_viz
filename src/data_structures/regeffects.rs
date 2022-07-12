use serde::{Deserialize, Serialize};

use crate::data_structures::DbID;
use rustc_hash::FxHashSet;

#[derive(Serialize, Deserialize, Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct Bucket(pub usize, pub u32);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RegEffectFacets(pub FxHashSet<DbID>, pub f32, pub f32);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RegEffectData {
    pub facets: Vec<RegEffectFacets>,
    pub associated_buckets: Vec<Bucket>,
    associated_buckets_set: FxHashSet<Bucket>,
}

impl RegEffectData {
    pub fn new() -> Self {
        let red = RegEffectData {
            facets: Vec::new(),
            associated_buckets: Vec::new(),
            associated_buckets_set: FxHashSet::default(),
        };
        red
    }

    pub fn add_facets(&mut self, facets: RegEffectFacets) {
        self.facets.push(facets);
    }

    pub fn update_buckets(&mut self, new_buckets: &FxHashSet<Bucket>) {
        for bucket in new_buckets - &self.associated_buckets_set {
            self.associated_buckets.push(bucket);
        }
        self.associated_buckets_set.extend(new_buckets.iter());
    }
}
