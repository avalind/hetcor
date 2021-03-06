
use clap::{Arg, App};
use rust_htslib::bcf::{IndexedReader, Read};
use std::fs::File;
use std::io::{self, BufRead};

#[derive(Debug)]
struct Region {
    ctg: String,
    start: u64,
    end: u64,
}

impl Region {
    fn from_string(line: &str) -> Result<Region, &str> {
        let mut parts = line.split("\t");
        
        let ctg = match parts.next() {
            Some(c) => c,
            None => return Err("Malformed record, while parsing contig"),
        };

        let start = match parts.next() {
            Some(s) => { 
                match s.parse::<u64>() {
                    Ok(i) => i,
                    Err(_e) => return Err("unable to parse to string"),   
                }
            },
            None => return Err("Malformed record while parsing start position"),
        };

        let end = match parts.next() {
            Some(s) => {
                match s.parse::<u64>() {
                    Ok(i) => i,
                    Err(_e) => return Err("unable to parse to string"),
                }
            },
            None => return Err("Malformed record while parsing end position"),
        };

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

fn main() {
    let matches = App::new("hetcor")
        .version("1.0")
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
        
        // first, find indices of the samples in question
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
        
        if regions == "" {
            for(_i, record_res) in bcf.records().enumerate() {
                let record = record_res.expect("malformed record!");
                let gts = record.genotypes().expect("Failed reading genotypes");
                let fgt = gts.get(fidx);
                let sgt = gts.get(sidx);
            
                if fgt == sgt && is_het(&format!("{}", fgt)) {
                    concordant_hets += 1;
                } 
            }
        } else {
            let bedfile = match File::open(regions) {
                Ok(file) => file,
                Err(_e) => {
                    panic!("Failed to open bed file!");
                },
            };
            
            let bedfile = io::BufReader::new(bedfile);
            let tally: i64 = bedfile.lines().map(|line| {                
                let unwrapped = line.unwrap();
                let reg = match Region::from_string(&unwrapped) {
                    Ok(r) => r,
                    Err(_e) => panic!("Malformed record passed to Region::from_string!"),
                };

                let hdr = bcf.header();
                let contig = match hdr.name2rid(reg.ctg.as_bytes()) {
                    Ok(cont) => cont,
                    Err(_e) => panic!("Unable to retrieve contig from bcf header!"),
                };
    
                // fetch the current region from the bcf file
                match bcf.fetch(contig, reg.start, reg.end) {
                    Ok(()) => (),
                    Err(_e) => panic!("Unable to fetch!"),
                }

                let subcount: i64 = bcf.records().map(|record| {
                    let site = match record {
                        Ok(s)  => s,
                        Err(_e) => panic!("Malformed vcf record"),
                    };

                    let gts = match site.genotypes() {
                        Ok(g) => g,
                        Err(_e) => panic!("genotype error!"),
                    };

                    if gts.get(fidx) == gts.get(sidx) && is_het(&format!("{}", gts.get(fidx))) {
                        return 1;
                    }
                    return 0;
                }).sum();
                return subcount;
            }).sum();
            concordant_hets = tally;
        }
        println!("{}", concordant_hets);

    }   
}   
