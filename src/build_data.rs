use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::rc::Rc;
use std::time::Instant;

use postgres::types::Json;
use postgres::{Client, Error};
use postgres_range::Range;

use crate::data_structures::facets::{
    facet_set, FACET_CCRE_CATEGORY, FACET_CCRE_OVERLAP, FACET_DIRECTION, FACET_EFFECT_SIZE,
    FACET_GRNA_TYPE, FACET_SIGNIFICANCE, FACET_TYPE_DISCRETE,
};
use crate::data_structures::*;
use crate::DbID;

const GRCH38: [(&str, i32); 25] = [
    ("1", 248956422),
    ("2", 242193529),
    ("3", 198295559),
    ("4", 190214555),
    ("5", 181538259),
    ("6", 170805979),
    ("7", 159345973),
    ("8", 145138636),
    ("9", 138394717),
    ("10", 133797422),
    ("11", 135086622),
    ("12", 133275309),
    ("13", 114364328),
    ("14", 107043718),
    ("15", 101991189),
    ("16", 90338345),
    ("17", 83257441),
    ("18", 80373285),
    ("19", 58617616),
    ("20", 64444167),
    ("21", 46709983),
    ("22", 50818468),
    ("X", 156040895),
    ("Y", 57227415),
    ("MT", 16569),
];
const GRCH37: [(&str, i32); 25] = [
    ("1", 249250621),
    ("2", 243199373),
    ("3", 198022430),
    ("4", 191154276),
    ("5", 180915260),
    ("6", 171115067),
    ("7", 159138663),
    ("8", 146364022),
    ("9", 141213431),
    ("10", 135534747),
    ("11", 135006516),
    ("12", 133851895),
    ("13", 115169878),
    ("14", 107349540),
    ("15", 102531392),
    ("16", 90354753),
    ("17", 81195210),
    ("18", 78077248),
    ("19", 59128983),
    ("20", 63025520),
    ("21", 48129895),
    ("22", 51304566),
    ("X", 155270560),
    ("Y", 59373566),
    ("MT", 16569),
];

fn select_assembly<'a>(
    assembly_name: &'a str,
    chromo: Option<&'a String>,
) -> Vec<(&'static str, i32)> {
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

    let assembly_info = select_assembly(options.assembly_name, options.chromo);

    let mut chrom_keys: HashMap<&str, usize> = HashMap::new();
    for i in 0..assembly_info.len() {
        chrom_keys.insert(assembly_info[i].0, i);
    }

    let mut source_buckets: HashMap<&str, Vec<Rc<RefCell<RegEffectData>>>> = HashMap::new();
    for chrom in &assembly_info {
        source_buckets.insert(
            chrom.0,
            (0..(bucket(chrom.1 as u32) + 1))
                .map(|_| Rc::new(RefCell::new(RegEffectData::new())))
                .collect(),
        );
    }
    let mut target_buckets: HashMap<&str, Vec<Rc<RefCell<RegEffectData>>>> = HashMap::new();
    for chrom in &assembly_info {
        target_buckets.insert(
            chrom.0,
            (0..(bucket(chrom.1 as u32) + 1))
                .map(|_| Rc::new(RefCell::new(RegEffectData::new())))
                .collect(),
        );
    }

    let mut chrom_data: Vec<ChromosomeData> = Vec::new();
    for chrom in &assembly_info {
        chrom_data.push(ChromosomeData::from(chrom.0, options.bucket_size))
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
        SELECT (search_regulatoryeffect_facet_values.regulatoryeffect_id) AS _prefetch_related_val_regulatoryeffect_id, search_facetvalue.id, search_facetvalue.value, search_facetvalue.facet_id
        FROM search_facetvalue
        INNER JOIN search_regulatoryeffect_facet_values ON (search_facetvalue.id = search_regulatoryeffect_facet_values.facetvalue_id)
        WHERE search_regulatoryeffect_facet_values.regulatoryeffect_id = ANY($1)"#
    )?;
    // (re id: DbID, dnaregion id: DbID, numeric facets: Json, chrom name: &str, location: Range(i32))
    let re_sources_statement = client.prepare(r#"
        SELECT (search_regulatoryeffect_sources.regulatoryeffect_id) AS _prefetch_related_val_regulatoryeffect_id, search_dnaregion.id, search_dnaregion.facet_num_values, search_dnaregion.chrom_name, search_dnaregion.location
        FROM search_dnaregion
        INNER JOIN search_regulatoryeffect_sources ON (search_dnaregion.id = search_regulatoryeffect_sources.dnaregion_id)
        WHERE search_regulatoryeffect_sources.regulatoryeffect_id = ANY($1)"#
    )?;
    // (re id: DbID, feature assembly id: DbID, chrom name: &str, location: Range(i32), strand: &str)
    let re_targets_statement = client.prepare(r#"
        SELECT (search_regulatoryeffect_target_assemblies.regulatoryeffect_id) AS _prefetch_related_val_regulatoryeffect_id, search_featureassembly.id, search_featureassembly.chrom_name, search_featureassembly.location, search_featureassembly.strand
        FROM search_featureassembly
        INNER JOIN search_regulatoryeffect_target_assemblies ON (search_featureassembly.id = search_regulatoryeffect_target_assemblies.featureassembly_id)
        WHERE search_regulatoryeffect_target_assemblies.regulatoryeffect_id = ANY($1)"#
    )?;
    // (source id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let source_facet_statement = client.prepare(r#"
        SELECT (search_dnaregion_facet_values.dnaregion_id) AS _prefetch_related_val_dnaregion_id, search_facetvalue.id, search_facetvalue.value, search_facetvalue.facet_id
        FROM search_facetvalue
        INNER JOIN search_dnaregion_facet_values ON (search_facetvalue.id = search_dnaregion_facet_values.facetvalue_id)
        WHERE search_dnaregion_facet_values.dnaregion_id = ANY($1)"#
    )?;
    let facet_range_statement = client.prepare(r#"
        SELECT MIN(((search_regulatoryeffect.facet_num_values -> $1))::double precision) AS min, MAX(((search_regulatoryeffect.facet_num_values -> $1))::double precision) AS max
        FROM search_regulatoryeffect
        INNER JOIN search_experiment ON (search_regulatoryeffect.experiment_id = search_experiment.id)
        WHERE search_experiment.accession_id = $2"#
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
    let source_facet_ids: HashSet<DbID> = HashSet::from([
        ccre_overlap_facet.id,
        ccre_category_facet.id,
        grna_type_facet.id,
    ]);

    let mut facet_ids: HashSet<DbID> = HashSet::new();

    let re_start_time = Instant::now();

    // (id: DbID, numeric facets: Json)
    let reg_effects_statement = client.prepare(r#"
        SELECT search_regulatoryeffect.id, search_regulatoryeffect.facet_num_values
        FROM search_regulatoryeffect
        INNER JOIN search_experiment ON (search_regulatoryeffect.experiment_id = search_experiment.id)
        WHERE search_experiment.accession_id = $1"#
    )?;
    let reg_effects_chromo_statement = client.prepare(r#"
        SELECT search_regulatoryeffect.id, search_regulatoryeffect.facet_num_values
        FROM search_regulatoryeffect
        INNER JOIN search_experiment ON (search_regulatoryeffect.experiment_id = search_experiment.id)
        INNER JOIN search_regulatoryeffect_target_assemblies as re_ta ON (search_regulatoryeffect.id = re_ta.regulatoryeffect_id)
        INNER JOIN search_featureassembly as fa ON (fa.id = re_ta.featureassembly_id)
        WHERE search_experiment.accession_id = $1 and fa.chrom_name = $2"#
    )?;
    let reg_effects = match options.chromo {
        None => client.query(&reg_effects_statement, &[&options.experiment_accession_id])?,
        Some(chromo) => client.query(
            &reg_effects_chromo_statement,
            &[&options.experiment_accession_id, &chromo],
        )?,
    };
    let mut reg_effect_num_facets: HashMap<DbID, HashMap<&str, f32>> = HashMap::new();
    for row in &reg_effects {
        let key = row.get::<usize, DbID>(0);
        let value = row.get::<usize, Json<HashMap<&str, f32>>>(1).0;
        reg_effect_num_facets.insert(key, value);
    }

    let reg_effect_id_list = reg_effects
        .iter()
        .map(|row| row.get::<&str, DbID>("id"))
        .collect::<Vec<DbID>>();
    let facet_values = client.query(&facet_values_statement, &[&reg_effect_id_list])?;
    // (re id: DbID, facet value id: DbID, value: &str, facet id: DbID)
    let mut facet_values_dict: HashMap<DbID, Vec<(DbID, &str, DbID)>> = HashMap::new();
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
    let mut source_dict: HashMap<
        DbID,
        Vec<(DbID, Option<Json<HashMap<&str, f32>>>, &str, Range<i32>)>,
    > = HashMap::new();
    for row in &sources {
        let key = row.get::<usize, DbID>(0);
        let value = (
            row.get::<usize, DbID>(1),
            row.get::<usize, Option<Json<HashMap<&str, f32>>>>(2),
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
    let mut target_dict: HashMap<DbID, Vec<(DbID, &str, Range<i32>, &str)>> = HashMap::new();
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
    let mut source_facet_dict: HashMap<DbID, Vec<(DbID, &str, DbID)>> = HashMap::new();
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

    let re_count: i64 = match options.chromo {
        None => client.query_one(r#"
            SELECT COUNT(*) AS count
            FROM search_regulatoryeffect
            INNER JOIN search_experiment ON (search_regulatoryeffect.experiment_id = search_experiment.id)
            WHERE search_experiment.accession_id = $1"#,
        &[&options.experiment_accession_id]
        )?,
        Some(chromo) => client.query_one(r#"
            SELECT COUNT(*) AS count
            FROM search_regulatoryeffect
            INNER JOIN search_experiment ON (search_regulatoryeffect.experiment_id = search_experiment.id)INNER JOIN search_regulatoryeffect_target_assemblies as re_ta ON (search_regulatoryeffect.id = re_ta.regulatoryeffect_id)
            INNER JOIN search_featureassembly as fa ON (fa.id = re_ta.featureassembly_id)
            WHERE search_experiment.accession_id = $1 and fa.chrom_name = $2"#,
        &[&options.experiment_accession_id, &chromo]
        )?,
    }.get("count");

    println!("Regulatory Effect count: {}", re_count);

    // For each regulatory effect we want to add all the facets associated with the effect itself,
    // its sources and its targets to the bucket associated with the each source and target.
    // For each source we want to keep track of all the target buckets it's associated with, and for each
    // source we want to keep track of all the source buckets it's associated with.
    for re_id in reg_effect_id_list {
        // If we are targeting a specific chromosome filter out sources not in that chromosome.
        // They won't be visible in the visualization and the bucket index will be wrong.
        let re_source = source_dict.get(&re_id).unwrap();
        let re_facets = reg_effect_num_facets.get(&re_id).unwrap();
        let effect_size = *re_facets.get(FACET_EFFECT_SIZE).unwrap();
        let significance = *re_facets.get(FACET_SIGNIFICANCE).unwrap();
        let all_chromos = match options.chromo {
            Some(_) => false,
            None => true,
        };
        let re_source_ref = Vec::from_iter(
            re_source
                .iter()
                .filter(|s| all_chromos || s.2 == options.chromo.unwrap()),
        );
        let mut source_counter: HashSet<Bucket> = HashSet::new();
        let mut target_counter: HashSet<Bucket> = HashSet::new();

        let mut source_disc_facets: HashSet<DbID> = HashSet::new();
        let mut reg_disc_facets: HashSet<DbID> = HashSet::new();

        for facets in facet_values_dict.get(&re_id) {
            facets
                .iter()
                .filter(|f| f.2 == dir_facet.id)
                .for_each(|f| drop(reg_disc_facets.insert(f.0)));
        }

        for source in &re_source_ref {
            source_counter.insert(Bucket(
                *chrom_keys
                    .get(source.2.strip_prefix("chr").unwrap())
                    .unwrap(),
                bucket(source.3.lower().unwrap().value as u32),
            ));

            for source_facets in &source_facet_dict.get(&source.0) {
                source_facets
                    .iter()
                    .filter(|f| source_facet_ids.contains(&f.2))
                    .for_each(|f| drop(source_disc_facets.insert(f.0)));
            }
        }

        let disc_facets = &reg_disc_facets | &source_disc_facets;

        for target in target_dict.get(&re_id).unwrap() {
            let chrom_name = target.1.strip_prefix("chr").unwrap();
            let target_start = match target.3 {
                "-" => target.2.upper().unwrap().value,
                _ => target.2.lower().unwrap().value,
            };
            let target_bucket = bucket(target_start as u32);
            target_counter.insert(Bucket(*chrom_keys.get(chrom_name).unwrap(), target_bucket));

            {
                let mut target_info = target_buckets
                    .get_mut(chrom_name)
                    .unwrap()
                    .get_mut(target_bucket as usize)
                    .unwrap()
                    .borrow_mut();
                target_info.add_facets(RegEffectFacets(
                    disc_facets.clone(),
                    effect_size,
                    significance,
                ));
                target_info.update_buckets(&source_counter);
            }
        }

        for source in &re_source_ref {
            let chrom_name = source.2.strip_prefix("chr").unwrap();

            let source_bucket = bucket(source.3.lower().unwrap().value as u32);

            {
                let mut source_info = source_buckets
                    .get_mut(chrom_name)
                    .unwrap()
                    .get_mut(source_bucket as usize)
                    .unwrap()
                    .borrow_mut();
                source_info.add_facets(RegEffectFacets(
                    disc_facets.clone(),
                    effect_size,
                    significance,
                ));
                source_info.update_buckets(&target_counter);
            }
        }

        facet_ids.extend(&disc_facets);
    }

    println!(
        "Buckets filled... {:.0}s",
        re_start_time.elapsed().as_secs()
    );

    for i in 0..assembly_info.len() {
        let chrom_name = assembly_info[i].0;
        let chrom_data = chrom_data.get_mut(i).unwrap();
        for j in 0..source_buckets.get(chrom_name).unwrap().len() {
            let source_bucket = source_buckets.get(chrom_name).unwrap().get(j).unwrap();
            if source_bucket.borrow().facets.len() == 0 {
                continue;
            }
            let source = Interval {
                start: options.bucket_size * (j as u32) + 1,
                values: Rc::clone(source_bucket),
            };
            chrom_data.source_intervals.push(source);
        }
        for j in 0..target_buckets.get(chrom_name).unwrap().len() {
            let target_bucket = target_buckets.get(chrom_name).unwrap().get(j).unwrap();
            if target_bucket.borrow().facets.len() == 0 {
                continue;
            }
            let target = Interval {
                start: options.bucket_size * (j as u32) + 1,
                values: Rc::clone(target_bucket),
            };
            chrom_data.target_intervals.push(target);
        }
    }

    // These are all the facets that are potentially relevant for coverage filtering
    let experiment_facet_coverages = facet_set();
    let experiment_facet_names: HashSet<&str> = HashSet::from([
        FACET_DIRECTION,
        FACET_EFFECT_SIZE,
        FACET_CCRE_CATEGORY,
        FACET_CCRE_OVERLAP,
        FACET_SIGNIFICANCE,
        FACET_GRNA_TYPE,
    ]);

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
        if facet.facet_type == FACET_TYPE_DISCRETE {
            let facet_values: HashMap<DbID, String> = all_facet_values
                .iter()
                .filter(|f| facet_ids.contains(&f.id))
                .map(|f| (f.id, f.value.to_string()))
                .collect();
            if facet_values.len() == 0 {
                continue;
            }
            facet.values = Some(facet_values);
        } else if facet.name == FACET_EFFECT_SIZE {
            let facet_range_row = client.query_one(
                &facet_range_statement,
                &[&FACET_EFFECT_SIZE, &options.experiment_accession_id],
            )?;
            facet.range = Some(FacetRange(
                facet_range_row.get("min"),
                facet_range_row.get("max"),
            ));
        } else if facet.name == FACET_SIGNIFICANCE {
            let facet_range_row = client.query_one(
                &facet_range_statement,
                &[&FACET_SIGNIFICANCE, &options.experiment_accession_id],
            )?;
            facet.range = Some(FacetRange(
                facet_range_row.get("min"),
                facet_range_row.get("max"),
            ));
        }

        facets.push(facet);
    }

    Ok(CoverageData {
        chromosomes: chrom_data,
        facets: facets.into_iter().map(|f| f.clone()).collect(),
    })
}
