use serde::{Deserialize, Serialize};

use rustc_hash::{FxHashMap, FxHashSet};

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
    pub coverage: Option<FxHashSet<FacetCoverage>>,
    pub range: Option<FacetRange>,
    pub values: Option<FxHashMap<DbID, String>>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct FacetValue {
    pub id: DbID,
    pub value: String,
    pub facet_id: DbID,
}

pub fn facet_set() -> FxHashMap<&'static str, FxHashSet<FacetCoverage>> {
    let mut experiment_facets: FxHashMap<&str, FxHashSet<FacetCoverage>> = FxHashMap::default();
    experiment_facets.insert(
        FACET_DIRECTION,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );
    experiment_facets.insert(
        FACET_EFFECT_SIZE,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );
    experiment_facets.insert(
        FACET_CCRE_CATEGORY,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );
    experiment_facets.insert(
        FACET_CCRE_OVERLAP,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );
    experiment_facets.insert(
        FACET_SIGNIFICANCE,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );
    experiment_facets.insert(
        FACET_GRNA_TYPE,
        FxHashSet::from_iter([FacetCoverage::Target, FacetCoverage::Source].into_iter()),
    );

    experiment_facets
}
