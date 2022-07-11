use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

use crate::data_structures::DbID;

pub const FACET_DIRECTION: &str = "Direction";
pub const FACET_EFFECT_SIZE: &str = "Effect Size";
pub const FACET_CCRE_CATEGORY: &str = "cCRE Category";
pub const FACET_CCRE_OVERLAP: &str = "cCRE Overlap";
pub const FACET_SIGNIFICANCE: &str = "Significance";
pub const FACET_GRNA_TYPE: &str = "gRNA Type";

pub const FACET_TYPE_DISCRETE: &str = "FacetType.DISCRETE";

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum FacetCoverage {
    Target,
    Source,
}

#[derive(Clone, Copy, Serialize, Deserialize, Debug)]
pub struct FacetRange(pub f32, pub f32);

#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct Facet {
    pub id: DbID,
    pub name: String,
    pub facet_type: String,
    pub description: String,
    pub coverage: Option<HashSet<FacetCoverage>>,
    pub range: Option<FacetRange>,
    pub values: Option<HashMap<DbID, String>>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FacetValue {
    pub id: DbID,
    pub value: String,
    pub facet_id: DbID,
}

pub fn facet_set() -> HashMap<&'static str, HashSet<FacetCoverage>> {
    let mut experiment_facets: HashMap<&str, HashSet<FacetCoverage>> = HashMap::new();
    experiment_facets.insert(
        FACET_DIRECTION,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );
    experiment_facets.insert(
        FACET_EFFECT_SIZE,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );
    experiment_facets.insert(
        FACET_CCRE_CATEGORY,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );
    experiment_facets.insert(
        FACET_CCRE_OVERLAP,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );
    experiment_facets.insert(
        FACET_SIGNIFICANCE,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );
    experiment_facets.insert(
        FACET_GRNA_TYPE,
        HashSet::from([FacetCoverage::Target, FacetCoverage::Source]),
    );

    experiment_facets
}
