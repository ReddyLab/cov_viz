use std::time::Instant;

use postgres::types::Json;
use postgres::{Client, Error};
use postgres_range::Range;
use roaring::RoaringTreemap;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::options::Options;

use cov_viz_ds::facets::{
    facet_set, FACET_CCRE_CATEGORY, FACET_CCRE_OVERLAP, FACET_DIRECTION, FACET_EFFECT_SIZE,
    FACET_GRNA_TYPE, FACET_SIGNIFICANCE, FACET_TYPE_CATEGORICAL,
};
use cov_viz_ds::*;

pub const MIN_SIG: f64 = 1e-100;

const GRCH38: [(&str, i32, u8); 25] = [
    ("1", 248956422, 0),
    ("2", 242193529, 1),
    ("3", 198295559, 2),
    ("4", 190214555, 3),
    ("5", 181538259, 4),
    ("6", 170805979, 5),
    ("7", 159345973, 6),
    ("8", 145138636, 7),
    ("9", 138394717, 8),
    ("10", 133797422, 9),
    ("11", 135086622, 10),
    ("12", 133275309, 11),
    ("13", 114364328, 12),
    ("14", 107043718, 13),
    ("15", 101991189, 14),
    ("16", 90338345, 15),
    ("17", 83257441, 16),
    ("18", 80373285, 17),
    ("19", 58617616, 18),
    ("20", 64444167, 19),
    ("21", 46709983, 20),
    ("22", 50818468, 21),
    ("X", 156040895, 22),
    ("Y", 57227415, 23),
    ("MT", 16569, 24),
];
const GRCH37: [(&str, i32, u8); 25] = [
    ("1", 249250621, 0),
    ("2", 243199373, 1),
    ("3", 198022430, 2),
    ("4", 191154276, 3),
    ("5", 180915260, 4),
    ("6", 171115067, 5),
    ("7", 159138663, 6),
    ("8", 146364022, 7),
    ("9", 141213431, 8),
    ("10", 135534747, 9),
    ("11", 135006516, 10),
    ("12", 133851895, 11),
    ("13", 115169878, 12),
    ("14", 107349540, 13),
    ("15", 102531392, 14),
    ("16", 90354753, 15),
    ("17", 81195210, 16),
    ("18", 78077248, 17),
    ("19", 59128983, 18),
    ("20", 63025520, 19),
    ("21", 48129895, 20),
    ("22", 51304566, 21),
    ("X", 155270560, 22),
    ("Y", 59373566, 23),
    ("MT", 16569, 24),
];

fn select_assembly<'a>(assembly_name: &'a str) -> Vec<(&'static str, i32, u8)> {
    match assembly_name {
        "GRCH37" => GRCH37.to_vec(),
        "GRCH38" => GRCH38.to_vec(),
        _ => {
            eprintln!(
                "Invalid genome {}. Must be \"GRCH37\" or \"GRCH38\"",
                assembly_name
            );
            panic!();
        }
    }
}

pub fn build_data(
    options: &Options,
    client: &mut Client,
) -> Result<(CoverageData, ExperimentFeatureData), Error> {
    let bucket = |size: u32| size / options.bucket_size;

    let assembly_info = select_assembly(&options.assembly_name);

    let mut significant_observations: Vec<ObservationData> = Vec::new();
    let mut nonsignificant_observations: Vec<ObservationData> = Vec::new();
    let mut feature_buckets = FxHashMap::<DbID, BucketLoc>::default();
    let mut source_set = RoaringTreemap::default();
    let mut target_set = RoaringTreemap::default();

    let mut chrom_keys: FxHashMap<&str, u8> = FxHashMap::default();
    for info in &assembly_info {
        chrom_keys.insert(info.0, info.2);
    }

    let chrom_data: Vec<ChromosomeData> = assembly_info
        .iter()
        .map(|chrom| ChromosomeData::from(chrom.0, chrom.2))
        .collect();

    let all_facet_rows = client.query(
        "SELECT id, name, description, facet_type FROM search_facet",
        &[],
    )?;
    let mut all_facets: Vec<Facet> = all_facet_rows
        .iter()
        .map(|r| Facet {
            id: r.get::<&str, i64>("id") as DbID,
            name: r.get::<&str, &str>("name").to_string(),
            description: r.get::<&str, &str>("description").to_string(),
            facet_type: r.get::<&str, &str>("facet_type").to_string(),
            coverage: None,
            range: None,
            range64: None,
            values: None,
        })
        .collect();

    let all_facet_value_rows =
        client.query("SELECT id, value, facet_id FROM search_facetvalue", &[])?;
    let all_facet_values: Vec<FacetValue> = all_facet_value_rows
        .iter()
        .map(|r| FacetValue {
            id: r.get::<&str, i64>("id") as DbID,
            value: r.get::<&str, &str>("value").to_string(),
            facet_id: r.get::<&str, i64>("facet_id") as DbID,
        })
        .collect();

    // (re id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let facet_values_statement = client.prepare(r#"
        SELECT (search_regulatoryeffectobservation_facet_values.regulatoryeffectobservation_id) AS _prefetch_related_val_regulatoryeffectobservation_id, search_facetvalue.id, search_facetvalue.value, search_facetvalue.facet_id
        FROM search_facetvalue
        INNER JOIN search_regulatoryeffectobservation_facet_values ON (search_facetvalue.id = search_regulatoryeffectobservation_facet_values.facetvalue_id)
        WHERE search_regulatoryeffectobservation_facet_values.regulatoryeffectobservation_id = ANY($1)"#
    )?;
    // (re id: DbID, dnafeature id: DbID, numeric facets: Json, chrom name: &str, location: Range(i32))
    let re_sources_statement = client.prepare(r#"
        SELECT (search_regulatoryeffectobservation_sources.regulatoryeffectobservation_id) AS _prefetch_related_val_regulatoryeffectobservation_id, search_dnafeature.id, search_dnafeature.facet_num_values, search_dnafeature.chrom_name, search_dnafeature.location
        FROM search_dnafeature
        INNER JOIN search_regulatoryeffectobservation_sources ON (search_dnafeature.id = search_regulatoryeffectobservation_sources.dnafeature_id)
        WHERE search_regulatoryeffectobservation_sources.regulatoryeffectobservation_id = ANY($1)"#
    )?;
    // (re id: DbID, feature assembly id: DbID, chrom name: &str, location: Range(i32), strand: &str)
    let re_targets_statement = client.prepare(r#"
        SELECT (search_regulatoryeffectobservation_targets.regulatoryeffectobservation_id) AS _prefetch_related_val_regulatoryeffectobservation_id, search_dnafeature.id, search_dnafeature.chrom_name, search_dnafeature.location, search_dnafeature.strand
        FROM search_dnafeature
        INNER JOIN search_regulatoryeffectobservation_targets ON (search_dnafeature.id = search_regulatoryeffectobservation_targets.dnafeature_id)
        WHERE search_regulatoryeffectobservation_targets.regulatoryeffectobservation_id = ANY($1)"#
    )?;
    // (source id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let source_facet_statement = client.prepare(r#"
        SELECT (search_dnafeature_facet_values.dnafeature_id) AS _prefetch_related_val_dnafeature_id, search_facetvalue.id, search_facetvalue.value, search_facetvalue.facet_id
        FROM search_facetvalue
        INNER JOIN search_dnafeature_facet_values ON (search_facetvalue.id = search_dnafeature_facet_values.facetvalue_id)
        WHERE search_dnafeature_facet_values.dnafeature_id = ANY($1)"#
    )?;
    let facet_range_statement = client.prepare(r#"
        SELECT MIN(((search_regulatoryeffectobservation.facet_num_values -> $1))::double precision) AS min, MAX(((search_regulatoryeffectobservation.facet_num_values -> $1))::double precision) AS max
        FROM search_regulatoryeffectobservation
        WHERE search_regulatoryeffectobservation.analysis_accession_id = $2"#
    )?;
    let dir_facet = all_facets
        .iter()
        .find(|f| f.name == FACET_DIRECTION)
        .unwrap();
    let ccre_overlap_facet = all_facets
        .iter()
        .find(|f| f.name == FACET_CCRE_OVERLAP)
        .unwrap();
    let ccre_category_facet = all_facets
        .iter()
        .find(|f| f.name == FACET_CCRE_CATEGORY)
        .unwrap();
    let grna_type_facet = all_facets
        .iter()
        .find(|f| f.name == FACET_GRNA_TYPE)
        .unwrap();
    let source_facet_ids: FxHashSet<DbID> = FxHashSet::from_iter(
        [
            ccre_overlap_facet.id,
            ccre_category_facet.id,
            grna_type_facet.id,
        ]
        .into_iter(),
    );

    let mut facet_ids: FxHashSet<DbID> = FxHashSet::default();

    let re_start_time = Instant::now();

    // (id: DbID, numeric facets: Json)
    let reg_effects_statement = client.prepare(r#"
        SELECT search_regulatoryeffectobservation.id, search_regulatoryeffectobservation.facet_num_values
        FROM search_regulatoryeffectobservation
        WHERE search_regulatoryeffectobservation.analysis_accession_id = $1"#
    )?;
    let reg_effects_chromo_statement = client.prepare(r#"
        SELECT search_regulatoryeffectobservation.id, search_regulatoryeffectobservation.facet_num_values
        FROM search_regulatoryeffectobservation
        INNER JOIN search_regulatoryeffectobservation_sources as re_s ON (search_regulatoryeffectobservation.id = re_s.regulatoryeffectobservation_id)
        INNER JOIN search_dnafeature as sf ON (sf.id = re_s.dnafeature_id)
        INNER JOIN search_regulatoryeffectobservation_targets as re_t ON (search_regulatoryeffectobservation.id = re_t.regulatoryeffectobservation_id)
        INNER JOIN search_dnafeature as tf ON (tf.id = re_t.dnafeature_id)
        WHERE search_regulatoryeffectobservation.analysis_accession_id = $1 and (sf.chrom_name = $2 or tf.chrom_name = $2)"#
    )?;
    let reg_effects = match &options.chromo {
        None => client.query(&reg_effects_statement, &[&options.analysis_accession_id])?,
        Some(chromo) => client.query(
            &reg_effects_chromo_statement,
            &[&options.analysis_accession_id, &chromo],
        )?,
    };
    let mut reg_effect_num_facets: FxHashMap<DbID, FxHashMap<&str, f32>> = FxHashMap::default();
    for row in &reg_effects {
        let key = row.get::<usize, i64>(0) as DbID;
        let value = row.get::<usize, Json<FxHashMap<&str, f32>>>(1).0;
        reg_effect_num_facets.insert(key, value);
    }

    let reg_effect_id_list = reg_effects
        .iter()
        .map(|row| row.get::<&str, i64>("id"))
        .collect::<Vec<i64>>();
    let facet_values = client.query(&facet_values_statement, &[&reg_effect_id_list])?;
    // (re id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let mut facet_values_dict: FxHashMap<DbID, Vec<(DbID, &str, DbID)>> = FxHashMap::default();
    for row in &facet_values {
        let key = row.get::<usize, i64>(0) as DbID;
        let value = (
            row.get::<usize, i64>(1) as DbID,
            row.get::<usize, &str>(2),
            row.get::<usize, i64>(3) as DbID,
        );
        facet_values_dict
            .entry(key)
            .and_modify(|e| e.push(value))
            .or_insert(vec![value]);
    }

    let sources = client.query(&re_sources_statement, &[&reg_effect_id_list])?;
    let mut source_dict: FxHashMap<
        DbID,
        Vec<(DbID, Option<Json<FxHashMap<&str, f32>>>, &str, Range<i32>)>,
    > = FxHashMap::default();
    for row in &sources {
        let key = row.get::<usize, i64>(0) as DbID;
        source_dict
            .entry(key)
            .and_modify(|e| {
                e.push((
                    row.get::<usize, i64>(1) as DbID,
                    row.get::<usize, Option<Json<FxHashMap<&str, f32>>>>(2),
                    row.get::<usize, &str>(3),
                    row.get::<usize, Range<i32>>(4),
                ))
            })
            .or_insert(vec![(
                row.get::<usize, i64>(1) as DbID,
                row.get::<usize, Option<Json<FxHashMap<&str, f32>>>>(2),
                row.get::<usize, &str>(3),
                row.get::<usize, Range<i32>>(4),
            )]);
    }

    let targets = client.query(&re_targets_statement, &[&reg_effect_id_list])?;
    let mut target_dict: FxHashMap<DbID, Vec<(DbID, &str, Range<i32>, &str)>> =
        FxHashMap::default();
    for row in &targets {
        let key = row.get::<usize, i64>(0) as DbID;
        let value = (
            row.get::<usize, i64>(1) as DbID,
            row.get::<usize, &str>(2),
            row.get::<usize, Range<i32>>(3),
            row.get::<usize, &str>(4),
        );
        target_dict
            .entry(key)
            .and_modify(|e| e.push(value))
            .or_insert(vec![value]);
    }

    let source_id_list = sources
        .iter()
        .map(|row| row.get::<&str, i64>("id"))
        .collect::<Vec<i64>>();
    let source_facets = client.query(&source_facet_statement, &[&source_id_list])?;
    let mut source_facet_dict: FxHashMap<DbID, Vec<(DbID, &str, DbID)>> = FxHashMap::default();
    for row in &source_facets {
        let key = row.get::<usize, i64>(0) as DbID;
        let value = (
            row.get::<usize, i64>(1) as DbID,
            row.get::<usize, &str>(2),
            row.get::<usize, i64>(3) as DbID,
        );
        source_facet_dict
            .entry(key)
            .and_modify(|e| e.push(value))
            .or_insert(vec![value]);
    }

    println!("Regulatory Effect count: {}", reg_effect_id_list.len());

    let nonsignificant_facet_value: DbID = all_facet_values
        .iter()
        .find(|fv| fv.facet_id == dir_facet.id && fv.value == "Non-significant")
        .unwrap()
        .id;

    // For each regulatory effect we want to add all the facets associated with the effect itself,
    // its sources and its targets to the bucket associated with the each source and target.
    // For each source we want to keep track of all the target buckets it's associated with, and for each
    // source we want to keep track of all the source buckets it's associated with.
    for reo_id in reg_effect_id_list {
        let re_facets = reg_effect_num_facets.get(&(reo_id as DbID)).unwrap();
        let effect_size = *re_facets.get(FACET_EFFECT_SIZE).unwrap();
        let significance: f64 = (*re_facets.get(FACET_SIGNIFICANCE).unwrap()).into();

        let re_sources = source_dict.get(&(reo_id as DbID)).unwrap();

        let mut source_counter: FxHashSet<BucketLoc> = FxHashSet::default();

        let mut source_cat_facets: FxHashSet<DbID> = FxHashSet::default();
        let mut reg_cat_facets: FxHashSet<DbID> = FxHashSet::default();

        // The only categorical REO facet we care about is the direction (depleted, enriched, or non-significant)
        if let Some(facets) = facet_values_dict.get(&(reo_id as DbID)) {
            facets
                .iter()
                .filter(|f| f.2 == dir_facet.id)
                .for_each(|f| drop(reg_cat_facets.insert(f.0)));
        }

        for source in re_sources {
            for source_facets in &source_facet_dict.get(&source.0) {
                source_facets
                    .iter()
                    .filter(|f| source_facet_ids.contains(&f.2))
                    .for_each(|f| drop(source_cat_facets.insert(f.0)));
            }

            let bucket_loc = BucketLoc {
                chrom: *chrom_keys
                    .get(source.2.strip_prefix("chr").unwrap())
                    .unwrap(),
                idx: bucket(source.3.lower().unwrap().value as u32),
            };
            source_counter.insert(bucket_loc);
            feature_buckets.insert(source.0, bucket_loc);
            source_set.insert(source.0);
        }

        let cat_facets = &reg_cat_facets | &source_cat_facets;

        let mut target_id: Option<DbID> = None;
        if let Some(targets) = target_dict.get(&(reo_id as DbID)) {
            let target = targets[0];
            target_id = Some(target.0);
            let chrom_name = target.1.strip_prefix("chr").unwrap();
            let x = chrom_keys.get(chrom_name);
            if let None = x {
                continue;
            }
            let target_start = match target.3 {
                "-" => target.2.upper().unwrap().value,
                _ => target.2.lower().unwrap().value,
            };
            let target_bucket = bucket(target_start as u32);
            let target_bucket = BucketLoc {
                chrom: *chrom_keys.get(chrom_name).unwrap(),
                idx: target_bucket,
            };
            feature_buckets.insert(target.0, target_bucket);
            target_set.insert(target.0);
        }

        if reg_cat_facets.contains(&nonsignificant_facet_value) {
            for (sid, _, _, _) in re_sources {
                nonsignificant_observations.push(ObservationData {
                    reo_id: reo_id as DbID,
                    facet_value_ids: cat_facets.iter().cloned().collect(),
                    source_id: *sid,
                    target_id,
                    effect_size,
                    significance,
                    neg_log_significance: -significance.max(MIN_SIG).log10(),
                });
            }
        } else {
            for (sid, _, _, _) in re_sources {
                significant_observations.push(ObservationData {
                    reo_id: reo_id as DbID,
                    facet_value_ids: cat_facets.iter().cloned().collect(),
                    source_id: *sid,
                    target_id,
                    effect_size,
                    significance,
                    neg_log_significance: -significance.max(MIN_SIG).log10(),
                });
            }
        }

        facet_ids.extend(&cat_facets);
    }

    println!(
        "Buckets filled... {:.0}s",
        re_start_time.elapsed().as_secs()
    );

    // These are all the facets that are potentially relevant for coverage filtering
    let experiment_facet_coverages = facet_set();
    let experiment_facet_names: FxHashSet<&str> = FxHashSet::from_iter(
        [
            FACET_DIRECTION,
            FACET_EFFECT_SIZE,
            FACET_CCRE_CATEGORY,
            FACET_CCRE_OVERLAP,
            FACET_SIGNIFICANCE,
            FACET_GRNA_TYPE,
        ]
        .into_iter(),
    );

    // The idea is to filter out facets that are in the database, but aren't used to annotate
    // data for this particular experiment.
    let mut facets = Vec::<&Facet>::new();
    for facet in all_facets
        .iter_mut()
        .filter(|f| experiment_facet_names.contains(f.name.as_str()))
    {
        facet.coverage = Some(
            experiment_facet_coverages
                .get(facet.name.as_str())
                .unwrap()
                .clone(),
        );
        if facet.facet_type == FACET_TYPE_CATEGORICAL {
            let facet_values: FxHashMap<DbID, String> = all_facet_values
                .iter()
                .filter(|f| facet_ids.contains(&f.id) && f.facet_id == facet.id)
                .map(|f| (f.id, f.value.to_string()))
                .collect();
            if facet_values.len() == 0 {
                continue;
            }
            facet.values = Some(facet_values);
        } else if facet.name == FACET_EFFECT_SIZE {
            let facet_range_row = client.query_one(
                &facet_range_statement,
                &[&FACET_EFFECT_SIZE, &options.analysis_accession_id],
            )?;
            facet.range = Some(FacetRange(
                facet_range_row.get::<&str, f64>("min") as f32,
                facet_range_row.get::<&str, f64>("max") as f32,
            ));
        } else if facet.name == FACET_SIGNIFICANCE {
            let facet_range_row = client.query_one(
                &facet_range_statement,
                &[&FACET_SIGNIFICANCE, &options.analysis_accession_id],
            )?;
            facet.range64 = Some(FacetRange64(
                facet_range_row.get::<&str, f64>("min"),
                facet_range_row.get::<&str, f64>("max"),
            ));
        }

        facets.push(facet);
    }

    Ok((
        CoverageData {
            significant_observations,
            nonsignificant_observations,
            bucket_size: options.bucket_size,
            chromosomes: chrom_data,
            facets: facets.into_iter().cloned().collect(),
            chrom_lengths: assembly_info.iter().map(|c| c.1 as usize).collect(),
            feature_buckets,
        },
        ExperimentFeatureData {
            sources: source_set,
            targets: target_set,
        },
    ))
}
