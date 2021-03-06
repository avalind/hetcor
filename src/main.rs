extern crate bio;

use std::env;
use clap::{Arg, App};


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
        println!("path to vcf file: {}", vcf_file);
        
        let first = match matches.value_of("first") {
            Some(val) => val, 
            None => panic!("No first sample specified!"),
        };

        let second = match matches.value_of("second") {
            Some(val) => val,
            None => panic!("No second sample specified!"),
        };
    
        
    }
    
    if let Some(regions) = matches.value_of("regions") {
        println!("path to region file: {}", regions);
    }
}
