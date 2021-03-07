
use clap::{Arg, App};
use rust_htslib::bcf::{IndexedReader, Read};
use anyhow::Result;
use std::fs::File;
use std::io::{self, BufRead};


#[derive(Debug)]
struct Region {
    ctg: String,
    start: u64,
    end: u64,
}

impl Region {
    fn from_string(line: &str) -> Result<Region> {
        let mut parts = line.split("\t");
        let ctg = parts.next().unwrap();
        let start = parts.next().unwrap().parse::<u64>()?;
        let end = parts.next().unwrap().parse::<u64>()?;
        return Ok(Region{
                ctg: ctg.to_string(),
                start: start,
                end: end})
    }
}

fn is_het(genotype: &str) -> bool {
    let accepted = ["1|0","0|1","1/0","0/1"];
    if accepted.iter().any(|&i| i == genotype) {
        return true;
    }
    return false;
}

fn run() -> Result<()> {
    let matches = App::new("hetcor")
        .version("0.1.0")
        .author("Anders Valind <anders.valind@med.lu.se>")
        .about("Calculates heterozygote concordance between two samples")
        .arg(Arg::new("regions")
            .short('r')
            .long("regions")
            .about("Specify the regions to process (.bed file)")
            .takes_value(true))
        .arg(Arg::new("vcf")
            .about("vcf file to process")
            .required(true)
            .index(3))
        .arg(Arg::new("first")
            .about("first sample in sample pair")
            .required(true)
            .index(1))
        .arg(Arg::new("second")
            .about("second sample in sample pair")
            .required(true)
            .index(2))
        .get_matches();

    if let Some(vcf_file) = matches.value_of("vcf") {
        let first = match matches.value_of("first") {
            Some(val) => val, 
            None => panic!("No first sample specified!"),
        };

        let second = match matches.value_of("second") {
            Some(val) => val,
            None => panic!("No second sample specified!"),
        };
        
        let regions = match matches.value_of("regions") {
            Some(r) => r,
            None => "",
        };
        
        let mut bcf = IndexedReader::from_path(vcf_file).expect("Unable to open vcf file");  
        let hdr = bcf.header();
        
        let fidx = match hdr.sample_id(first.as_bytes()) {
            Some(idx) => idx,
            None => {
                panic!("Sample {} not present in the vcf file header", first);
            },
        };

        let sidx = match hdr.sample_id(second.as_bytes()) {
            Some(idx) => idx,
            None => {
                panic!("Sample {} not present in the vcf file header", second);
            },
        };
        
        let mut concordant_hets = 0;
        let mut total_variants = 0;
        if regions == "" {
            for(_i, record_res) in bcf.records().enumerate() {
                let record = record_res?;
                let gts = record.genotypes()?;
                let fgt = gts.get(fidx);
                let sgt = gts.get(sidx);
            
                if fgt == sgt && is_het(&format!("{}", fgt)) {
                    concordant_hets += 1;
                }; 
                total_variants += 1;
            }
        } else {
            let bedfile = File::open(regions)?;
            let bedfile = io::BufReader::new(bedfile);
            let tally: Result<i64> = bedfile.lines().map(|line| {                
                let unwrapped = line.unwrap();
                let reg = Region::from_string(&unwrapped)?;
                let hdr = bcf.header();
                let contig = hdr.name2rid(reg.ctg.as_bytes())?;
    
                // fetch the current region from the bcf file
                bcf.fetch(contig, reg.start, reg.end)?;
                let subcount: Result<i64> = bcf.records().map(|record| {
                    let site = record?;
                    let gts = site.genotypes()?;
                    total_variants += 1;
                    if gts.get(fidx) == gts.get(sidx) && is_het(&format!("{}", gts.get(fidx))) {
                        return Ok(1)
                    }
                    return Ok(0);
                }).sum();
                
                return subcount;
            }).sum();
            concordant_hets = tally?;
        }
        println!("total_variants\tconcordant_hets");
        println!("{}\t{}", total_variants, concordant_hets);
    }   
    return Ok(());
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error while executing: {}", e);
        e.chain().skip(1).for_each(|cause| eprintln!("because: {}", cause));
        std::process::exit(1);
    }
}   
