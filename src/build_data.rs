use std::time::Instant;

use postgres::types::Json;
use postgres::{Client, Error};
use postgres_range::Range;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::options::Options;

use cov_viz_ds::facets::{
    facet_set, FACET_CCRE_CATEGORY, FACET_CCRE_OVERLAP, FACET_DIRECTION, FACET_EFFECT_SIZE,
    FACET_GRNA_TYPE, FACET_SIGNIFICANCE, FACET_TYPE_CATEGORICAL,
};
use cov_viz_ds::*;

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

fn select_assembly<'a>(
    assembly_name: &'a str,
    chromo: &Option<String>,
) -> Vec<(&'static str, i32, u8)> {
    let mut assembly_info = match assembly_name {
        "GRCH37" => GRCH37.to_vec(),
        "GRCH38" => GRCH38.to_vec(),
        _ => {
            eprintln!(
                "Invalid genome {}. Must be \"GRCH37\" or \"GRCH38\"",
                assembly_name
            );
            panic!();
        }
    };

    if let Some(chromo) = chromo {
        let chromo = chromo.strip_prefix("chr").unwrap();
        assembly_info = assembly_info
            .into_iter()
            .filter(|c| c.0 == chromo)
            .collect();
    }

    assembly_info
}

pub fn build_data(options: &Options, client: &mut Client) -> Result<CoverageData, Error> {
    let bucket = |size: u32| size / options.bucket_size;

    let assembly_info = select_assembly(&options.assembly_name, &options.chromo);

    let mut chrom_keys: FxHashMap<&str, u8> = FxHashMap::default();
    for i in 0..assembly_info.len() {
        chrom_keys.insert(assembly_info[i].0, assembly_info[i].2);
    }

    let mut source_buckets: FxHashMap<&str, Vec<FxHashMap<DbID, FeatureData>>> =
        FxHashMap::default();
    for chrom in &assembly_info {
        source_buckets.insert(
            chrom.0,
            (0..(bucket(chrom.1 as u32) + 1))
                .map(|_| FxHashMap::default())
                .collect(),
        );
    }
    let mut target_buckets: FxHashMap<&str, Vec<FxHashMap<DbID, FeatureData>>> =
        FxHashMap::default();
    for chrom in &assembly_info {
        target_buckets.insert(
            chrom.0,
            (0..(bucket(chrom.1 as u32) + 1))
                .map(|_| FxHashMap::default())
                .collect(),
        );
    }

    let mut chrom_data: Vec<ChromosomeData> = Vec::new();
    for chrom in &assembly_info {
        chrom_data.push(ChromosomeData::from(chrom.0, chrom.2, options.bucket_size))
    }

    let all_facet_rows = client.query(
        "SELECT id, name, description, facet_type FROM search_facet",
        &[],
    )?;
    let mut all_facets: Vec<Facet> = all_facet_rows
        .iter()
        .map(|r| Facet {
            id: r.get::<&str, DbID>("id"),
            name: r.get::<&str, &str>("name").to_string(),
            description: r.get::<&str, &str>("description").to_string(),
            facet_type: r.get::<&str, &str>("facet_type").to_string(),
            coverage: None,
            range: None,
            values: None,
        })
        .collect();

    let all_facet_value_rows =
        client.query("SELECT id, value, facet_id FROM search_facetvalue", &[])?;
    let all_facet_values: Vec<FacetValue> = all_facet_value_rows
        .iter()
        .map(|r| FacetValue {
            id: r.get::<&str, DbID>("id"),
            value: r.get::<&str, &str>("value").to_string(),
            facet_id: r.get::<&str, DbID>("facet_id"),
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
        INNER JOIN search_regulatoryeffectobservation_sources as re_ta ON (search_regulatoryeffectobservation.id = re_ta.regulatoryeffectobservation_id)
        INNER JOIN search_dnafeature as fa ON (fa.id = re_ta.dnafeature_id)
        WHERE search_regulatoryeffectobservation.analysis_accession_id = $1 and fa.chrom_name = $2"#
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
        let key = row.get::<usize, DbID>(0);
        let value = row.get::<usize, Json<FxHashMap<&str, f32>>>(1).0;
        reg_effect_num_facets.insert(key, value);
    }

    let reg_effect_id_list = reg_effects
        .iter()
        .map(|row| row.get::<&str, DbID>("id"))
        .collect::<Vec<DbID>>();
    let facet_values = client.query(&facet_values_statement, &[&reg_effect_id_list])?;
    // (re id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let mut facet_values_dict: FxHashMap<DbID, Vec<(DbID, &str, DbID)>> = FxHashMap::default();
    for row in &facet_values {
        let key = row.get::<usize, DbID>(0);
        let value = (
            row.get::<usize, DbID>(1),
            row.get::<usize, &str>(2),
            row.get::<usize, DbID>(3),
        );
        if facet_values_dict.contains_key(&key) {
            if let Some(v) = facet_values_dict.get_mut(&key) {
                v.push(value);
            }
        } else {
            facet_values_dict.insert(key, vec![value]);
        }
    }

    let sources = client.query(&re_sources_statement, &[&reg_effect_id_list])?;
    let mut source_dict: FxHashMap<
        DbID,
        Vec<(DbID, Option<Json<FxHashMap<&str, f32>>>, &str, Range<i32>)>,
    > = FxHashMap::default();
    for row in &sources {
        let key = row.get::<usize, DbID>(0);
        let value = (
            row.get::<usize, DbID>(1),
            row.get::<usize, Option<Json<FxHashMap<&str, f32>>>>(2),
            row.get::<usize, &str>(3),
            row.get::<usize, Range<i32>>(4),
        );
        if source_dict.contains_key(&key) {
            if let Some(v) = source_dict.get_mut(&key) {
                v.push(value);
            }
        } else {
            source_dict.insert(key, vec![value]);
        }
    }

    let targets = client.query(&re_targets_statement, &[&reg_effect_id_list])?;
    let mut target_dict: FxHashMap<DbID, Vec<(DbID, &str, Range<i32>, &str)>> =
        FxHashMap::default();
    for row in &targets {
        let key = row.get::<usize, DbID>(0);
        let value = (
            row.get::<usize, DbID>(1),
            row.get::<usize, &str>(2),
            row.get::<usize, Range<i32>>(3),
            row.get::<usize, &str>(4),
        );
        if target_dict.contains_key(&key) {
            if let Some(v) = target_dict.get_mut(&key) {
                v.push(value);
            }
        } else {
            target_dict.insert(key, vec![value]);
        }
    }

    let source_id_list = sources
        .iter()
        .map(|row| row.get::<&str, DbID>("id"))
        .collect::<Vec<DbID>>();
    let source_facets = client.query(&source_facet_statement, &[&source_id_list])?;
    let mut source_facet_dict: FxHashMap<DbID, Vec<(DbID, &str, DbID)>> = FxHashMap::default();
    for row in &source_facets {
        let key = row.get::<usize, DbID>(0);
        let value = (
            row.get::<usize, DbID>(1),
            row.get::<usize, &str>(2),
            row.get::<usize, DbID>(3),
        );
        if source_facet_dict.contains_key(&key) {
            if let Some(v) = source_facet_dict.get_mut(&key) {
                v.push(value);
            }
        } else {
            source_facet_dict.insert(key, vec![value]);
        }
    }

    println!("Regulatory Effect count: {}", reg_effect_id_list.len());

    // For each regulatory effect we want to add all the facets associated with the effect itself,
    // its sources and its targets to the bucket associated with the each source and target.
    // For each source we want to keep track of all the target buckets it's associated with, and for each
    // source we want to keep track of all the source buckets it's associated with.
    for reo_id in reg_effect_id_list {
        // If we are targeting a specific chromosome filter out sources not in that chromosome.
        // They won't be visible in the visualization and the bucket index will be wrong.
        let re_source = source_dict.get(&reo_id).unwrap();
        let re_facets = reg_effect_num_facets.get(&reo_id).unwrap();
        let effect_size = *re_facets.get(FACET_EFFECT_SIZE).unwrap();
        let significance = *re_facets.get(FACET_SIGNIFICANCE).unwrap();
        let all_chromos = match options.chromo {
            Some(_) => false,
            None => true,
        };
        let re_source_ref = Vec::from_iter(
            re_source
                .iter()
                .filter(|s| all_chromos || s.2 == options.chromo.as_ref().unwrap()),
        );
        let mut source_counter: FxHashSet<BucketLoc> = FxHashSet::default();
        let mut target_counter: FxHashSet<BucketLoc> = FxHashSet::default();

        let mut source_disc_facets: FxHashSet<DbID> = FxHashSet::default();
        let mut reg_disc_facets: FxHashSet<DbID> = FxHashSet::default();

        if let Some(facets) = facet_values_dict.get(&reo_id) {
            facets
                .iter()
                .filter(|f| f.2 == dir_facet.id)
                .for_each(|f| drop(reg_disc_facets.insert(f.0)));
        }

        for source in &re_source_ref {
            source_counter.insert(BucketLoc {
                chrom: *chrom_keys
                    .get(source.2.strip_prefix("chr").unwrap())
                    .unwrap(),
                idx: bucket(source.3.lower().unwrap().value as u32),
            });

            for source_facets in &source_facet_dict.get(&source.0) {
                source_facets
                    .iter()
                    .filter(|f| source_facet_ids.contains(&f.2))
                    .for_each(|f| drop(source_disc_facets.insert(f.0)));
            }
        }

        let disc_facets = &reg_disc_facets | &source_disc_facets;

        if let Some(targets) = target_dict.get(&reo_id) {
            for target in targets {
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
                target_counter.insert(BucketLoc {
                    chrom: *chrom_keys.get(chrom_name).unwrap(),
                    idx: target_bucket,
                });

                {
                    let targets = target_buckets
                        .get_mut(chrom_name)
                        .unwrap()
                        .get_mut(target_bucket as usize)
                        .unwrap();
                    let target_info = targets
                        .entry(target.0)
                        .or_insert(FeatureData::new(target.0));
                    target_info.add_facets(ObservationData {
                        reo_id,
                        facet_ids: disc_facets.clone().into_iter().collect(),
                        effect_size,
                        significance,
                    });
                    target_info.update_buckets(&source_counter);
                }
            }
        }

        for source in &re_source_ref {
            let chrom_name = source.2.strip_prefix("chr").unwrap();

            let source_bucket = bucket(source.3.lower().unwrap().value as u32);

            {
                let sources = source_buckets
                    .get_mut(chrom_name)
                    .unwrap()
                    .get_mut(source_bucket as usize)
                    .unwrap();
                let source_info = sources
                    .entry(source.0)
                    .or_insert(FeatureData::new(source.0));
                source_info.add_facets(ObservationData {
                    reo_id,
                    facet_ids: disc_facets.clone().into_iter().collect(),
                    effect_size,
                    significance,
                });
                source_info.update_buckets(&target_counter);
            }
        }

        facet_ids.extend(&disc_facets);
    }

    println!(
        "Buckets filled... {:.0}s",
        re_start_time.elapsed().as_secs()
    );

    for (i, info) in assembly_info.iter().enumerate() {
        let chrom_name = info.0;
        let chrom_data = chrom_data.get_mut(i).unwrap();
        for (j, source_bucket) in source_buckets.get(chrom_name).unwrap().iter().enumerate() {
            {
                if source_bucket.len() == 0 {
                    continue;
                }
            }
            let source = Interval {
                start: options.bucket_size * (j as u32) + 1,
                values: source_bucket.clone().into_values().collect(),
            };
            chrom_data.source_intervals.push(source);
        }
        for (j, target_bucket) in target_buckets.get(chrom_name).unwrap().iter().enumerate() {
            {
                if target_bucket.len() == 0 {
                    continue;
                }
            }
            let target = Interval {
                start: options.bucket_size * (j as u32) + 1,
                values: target_bucket.clone().into_values().collect(),
            };
            chrom_data.target_intervals.push(target);
        }
    }

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
            facet.range = Some(FacetRange(
                facet_range_row.get::<&str, f64>("min") as f32,
                facet_range_row.get::<&str, f64>("max") as f32,
            ));
        }

        facets.push(facet);
    }

    Ok(CoverageData {
        chromosomes: chrom_data,
        facets: facets.into_iter().map(|f| f.clone()).collect(),
        chrom_lengths: assembly_info.iter().map(|c| c.1 as usize).collect(),
    })
}
