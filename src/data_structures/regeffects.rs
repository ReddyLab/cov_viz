use std::collections::HashSet;
use serde::{Deserialize, Serialize};

use crate::DbID;

#[derive(Serialize, Deserialize, Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub struct Bucket(pub usize, pub u32);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RegEffectFacets(pub HashSet<DbID>, pub f32, pub f32);

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct RegEffectData {
    pub facets: Vec<RegEffectFacets>,
    pub associated_buckets: HashSet<Bucket>,
}

impl RegEffectData {
    pub fn new() -> Self {
        let red = RegEffectData {
            facets: Vec::new(),
            associated_buckets: HashSet::new(),
        };
        red
    }

    pub fn add_facets(&mut self, facets: RegEffectFacets) {
        self.facets.push(facets);
    }

    pub fn update_buckets(&mut self, new_buckets: &HashSet<Bucket>) {
        for bucket in new_buckets {
            self.associated_buckets.insert(*bucket);
        }
    }
}
