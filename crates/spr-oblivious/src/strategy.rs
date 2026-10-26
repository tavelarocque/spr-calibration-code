//! An oblivious no-reuse strategy: a finite list of tuples of distinct cells with
//! exact rational probabilities. Serialised as JSON.

use crate::rational::{fmt, parse};
use num_rational::BigRational;
use serde_json::{Value, json};

#[derive(Clone, Debug)]
pub struct Strategy {
    pub n: usize,
    pub s: usize,
    /// `(tuple, probability)`; cells are 0-based.
    pub entries: Vec<(Vec<u8>, BigRational)>,
}

impl Strategy {
    pub fn to_json(&self, extra: Value) -> Value {
        let tuples: Vec<Value> =
            self.entries.iter().map(|(t, p)| json!([t.iter().map(|&c| c as u64).collect::<Vec<_>>(), fmt(p)])).collect();
        json!({
            "n": self.n,
            "s": self.s,
            "cells": "0-based; a tuple lists the cells in the order the Pointer visits them",
            "info": extra,
            "tuples": tuples,
        })
    }

    pub fn from_json(v: &Value) -> Result<Self, String> {
        let n = v["n"].as_u64().ok_or("missing n")? as usize;
        let s = v["s"].as_u64().ok_or("missing s")? as usize;
        let mut entries = Vec::new();
        for e in v["tuples"].as_array().ok_or("missing tuples")? {
            let t: Vec<u8> = e[0]
                .as_array()
                .ok_or("bad tuple")?
                .iter()
                .map(|c| c.as_u64().map(|c| c as u8).ok_or("bad cell"))
                .collect::<Result<_, _>>()?;
            let p = parse(e[1].as_str().ok_or("bad probability")?).ok_or("bad rational")?;
            entries.push((t, p));
        }
        Ok(Strategy { n, s, entries })
    }
}
