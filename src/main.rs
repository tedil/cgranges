use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::ops::{Deref, Range};
use std::path::PathBuf;

use bio::data_structures::interval_tree::ArrayBackedIntervalTree;
use csv;
use serde_derive::{Deserialize, Serialize};
use structopt::StructOpt;

#[derive(StructOpt)]
#[structopt(name = "bedcov")]
struct Opt {
    #[structopt(parse(from_os_str))]
    a: PathBuf,

    #[structopt(parse(from_os_str))]
    b: PathBuf,
}

#[derive(Debug, Serialize, Default, Deserialize, Clone)]
pub struct Record {
    chrom: String,
    start: u64,
    end: u64,
}

fn main() -> Result<(), std::io::Error> {
    let args = Opt::from_args();
    let bed_ref_reader = BufReader::new(File::open(args.a)?);
    let mut bed_ref = csv::ReaderBuilder::default()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(bed_ref_reader);

    let mut trees: HashMap<String, ArrayBackedIntervalTree<_, _>> = HashMap::with_capacity(24);
    bed_ref
        .deserialize()
        .filter_map(Result::ok)
        .for_each(|r: Record| {
            trees
                .entry(r.chrom)
                .or_insert(ArrayBackedIntervalTree::new())
                .insert(r.start..r.end, 0)
        });
    // sorting happens in `index()`
    trees.values_mut().for_each(|tree| tree.index());

    let bed_target_reader = BufReader::new(File::open(args.b)?);
    let mut bed_target = csv::ReaderBuilder::default()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(bed_target_reader);
    let mut buffer = Vec::with_capacity(4096);
    bed_target
        .deserialize()
        .filter_map(Result::ok)
        .for_each(|r: Record| {
            let (target_start, target_end, chrom) = (r.start, r.end, r.chrom.to_string());
            let (mut cov, mut cov_start, mut cov_end, mut n) = (0, 0, 0, 0);
            trees[&chrom].find_into(target_start..target_end, &mut buffer);
            for entry in &buffer {
                n += 1;
                let Range { mut start, mut end } = entry.interval().deref();
                start = start.max(target_start);
                end = end.min(target_end);
                if start > cov_end {
                    cov += cov_end - cov_start;
                    cov_start = start;
                    cov_end = end;
                } else if cov_end < end {
                    cov_end = end;
                }
            }
            cov += cov_end - cov_start;
            println!(
                "{chrom}\t{start}\t{end}\t{count}\t{coverage}",
                chrom = chrom,
                start = target_start,
                end = target_end,
                count = n,
                coverage = cov
            );
        });
    Ok(())
}
