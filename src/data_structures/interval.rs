use std::collections::HashMap;
use std::fmt;
use std::ops::Deref;
use std::sync::{Arc, Mutex};

use serde::de::{self, Deserializer, MapAccess, SeqAccess, Visitor};
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};

use crate::data_structures::RegEffectData;
use crate::DbID;

#[derive(Debug)]
pub struct Interval {
    pub start: u32,
    pub values: Arc<Mutex<HashMap<DbID, RegEffectData>>>,
}

const INTERVAL_FIELD_START: &str = "start";
const INTERVAL_FIELD_VALUES: &str = "values";

impl Interval {
    fn new(start: u32, values: HashMap<DbID, RegEffectData>) -> Self {
        let i = Interval {
            start: start,
            values: Arc::new(Mutex::new(values)),
        };

        i
    }
}

impl Serialize for Interval {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Interval", 2)?;
        state.serialize_field(INTERVAL_FIELD_START, &self.start)?;
        {
            let values = self.values.lock().unwrap();
            state.serialize_field(INTERVAL_FIELD_VALUES, values.deref())?;
        }

        state.end()
    }
}

impl<'de> Deserialize<'de> for Interval {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            Start,
            Values,
        }

        struct IntervalVisitor;

        impl<'de> Visitor<'de> for IntervalVisitor {
            type Value = Interval;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct Interval")
            }

            fn visit_seq<V>(self, mut seq: V) -> Result<Interval, V::Error>
            where
                V: SeqAccess<'de>,
            {
                let start = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(0, &self))?;
                let values = seq
                    .next_element()?
                    .ok_or_else(|| de::Error::invalid_length(1, &self))?;
                Ok(Interval::new(start, values))
            }

            fn visit_map<V>(self, mut map: V) -> Result<Interval, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut start = None;
                let mut values = None;
                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Start => {
                            if start.is_some() {
                                return Err(de::Error::duplicate_field(INTERVAL_FIELD_START));
                            }
                            start = Some(map.next_value()?);
                        }
                        Field::Values => {
                            if values.is_some() {
                                return Err(de::Error::duplicate_field(INTERVAL_FIELD_VALUES));
                            }
                            values = Some(map.next_value()?);
                        }
                    }
                }
                let start = start.ok_or_else(|| de::Error::missing_field(INTERVAL_FIELD_START))?;
                let values =
                    values.ok_or_else(|| de::Error::missing_field(INTERVAL_FIELD_VALUES))?;
                Ok(Interval::new(start, values))
            }
        }

        const FIELDS: &'static [&'static str] = &[INTERVAL_FIELD_START, INTERVAL_FIELD_VALUES];
        deserializer.deserialize_struct("Interval", FIELDS, IntervalVisitor)
    }
}
