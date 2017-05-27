#![feature(test)]

extern crate test;
extern crate fastx;
extern crate bio;

use std::io::BufReader;
use std::fs::File;
use std::io::{Cursor, Read};

use fastx::Records;
use bio::io::{fasta, fastq};

#[bench]
fn fastx_fasta(b: &mut test::Bencher) {
    let mut data = String::new();
    let mut file = File::open("tests/test_data/test.fa").unwrap();
    file.read_to_string(&mut data);
        
    b.iter(|| {
        let records = Records::from_fasta(data.as_bytes());
        let mut x = 0;
        for record in records {
            x += record.seq.len();
        }
    })
}

#[bench]
fn rustbio_fasta(b: &mut test::Bencher) {
    let mut data = String::new();
    let mut file = File::open("tests/test_data/test.fa").unwrap();
    file.read_to_string(&mut data);
    b.iter(|| {
        let reader = fasta::Reader::new(data.as_bytes());
        let mut x = 0;
        for record in reader.records() {
            x += record.unwrap().seq().len();
        }
    })
}


#[bench]
fn fastx_fastq(b: &mut test::Bencher) {
    let mut data = String::new();
    let mut file = File::open("tests/test_data/test.fq").unwrap();
    file.read_to_string(&mut data);
    b.iter(|| {
        let records = Records::from_fasta(data.as_bytes());
        let mut x = 0;
        for record in records {
            x += record.seq.len();
        }
    })
}

#[bench]
fn rustbio_fastq(b: &mut test::Bencher) {
    let mut data = String::new();
    let mut file = File::open("tests/test_data/test.fq").unwrap();
    file.read_to_string(&mut data);
    b.iter(|| {
        let reader = fastq::Reader::new(data.as_bytes());
        let mut x = 0;
        for record in reader.records() {
            x += record.unwrap().seq().len();
        }
    })
}