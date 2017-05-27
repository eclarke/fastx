#[macro_use]
extern crate nom;
extern crate itertools;

use std::io;

use nom::IResult;

mod internals;
use internals::*;
pub use internals::Record;

const FA_REC: u8 = '>' as u8;

pub struct Records<R: io::BufRead> {
    input: R,
    parser: fn(&mut R) -> Option<Record>
}

impl<'a, R: io::BufRead> Records<R> {
    pub fn from_fasta(input: R) -> Records<R> {
        let mut input = input;
        input.read_until(FA_REC, &mut Vec::new()).expect("fasta parser: error while reading input!");
        Records {
            input, 
            parser: read_fasta
        }
    }

    pub fn from_fastq(input: R) -> Records<R> {
        Records {
            input,
            parser: read_fastq
        }
    }
}

impl<R: io::BufRead> Iterator for Records<R> {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        (self.parser)(&mut self.input)
    }
}

fn read_fasta<R: io::BufRead>(reader: &mut R) -> Option<Record> {
    let mut buf = Vec::new();
    reader.read_until(FA_REC, &mut buf).expect("fasta parser: error while reading from input!");
    buf.insert(0, FA_REC);

    match parse_single_fasta(&buf[..buf.len()]) {
        IResult::Done(_, record) => Some(record),
        _ => None,
    }
}

// #[allow(dead_code)]
fn read_fastq<R: io::BufRead>(reader: &mut R) -> Option<Record> {
    let mut buf: Vec<u8> = Vec::new();
    loop {
        let ibuf_len: usize;
        {
            let ibuf = reader.fill_buf().expect("fastq parser: error while reading from fastq file!");  
            ibuf_len = ibuf.len();
            if ibuf_len == 0 {
                return None
            }
            buf.extend(ibuf);
        }
        match parse_single_fastq(&buf) {
            IResult::Done(extra, record) => {
                // Consume only as much as was needed to get the record
                reader.consume(ibuf_len - extra.len());
                return Some(record)
            },
            IResult::Error(_) => return None,
            IResult::Incomplete(_) => {
                reader.consume(ibuf_len);
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::BufReader;
    use super::*;
    #[test]
    fn test_fasta_reader() {
        let fa = BufReader::new(File::open("tests/test_data/test.fa").unwrap());
        let mut records = Records::from_fasta(fa);
        let rec1 = records.next();
        let rec2 = records.next();
        let rec3 = records.next();
        assert_eq!(None, records.next());
        assert_eq!(
            rec1, 
            Some(Record{
                id: String::from("gi|1348912|gb|G26680|G26680 human STS STS_D11729.^Agi|1396336|gb|G27617|G27617 human STS SHGC-32648."),
                seq: b"CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANATNATGACTCTACTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGAGAAACAAACAAGCCCCACTCAGCTCAGTCACAGGAGAGANCACAGAAAGTCTTAGGATCATGANCTCTGAAAAAAAGAGAAACCTTATCTTTNCTTTGTGGTTCCTTTAAACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATGTNAAGCNTGACAACACCAGGCAGGTATGAGAGGAAAGCAAGGCCCGTCCATNGCTGTCCGTACNCTTACGGNTTGCTTGTNGGAGNCATTTNGGTATTGTTTGTTGTAANANCCAAAANGGGCTTTGGNNTGGNAAAAGGGCAGANNGGGGGGGTTGGTGTNGTTTTTTGGGGGGANNNTTTNGATTTGGTNCCGGGNTTTNGTTTNCCNCGGNACCGGNTTTTGGTTGGGGNCCATTTNTGNGGGGCNTTGGNGTTNCNTTNCCCNNNTNNGANTGGTTTNA"
                    .to_vec(),
                qual: None
            })
        );
        assert_eq!(
            rec2.unwrap().id, String::from("gi|1348917|gb|G26685|G26685 human STS STS_D11734.")
        );
        assert_eq!(
            rec3.unwrap().seq.len(), 471
        );
    }

    #[test]
    fn test_fastq_reader() {
        let fq = BufReader::new(File::open("tests/test_data/test.fq").unwrap());
        let mut records = Records::from_fastq(fq);
        let rec1 = records.next();
        let rec2 = records.next();
        let rec3 = records.next();
        assert_eq!(None, records.next());
        assert_eq!(
            rec1,
            Some(Record{
                id: String::from("D00727:24:C90LYANXX:7:1109:1635:2137 1:N:0:AAGAGGCA+GTAAGGAG"),
                seq: b"GTCCAGAGCTTCGGTATAACGCTTGATCGCCAATCATTTTCGGCGCAGGATCACTCGATGAGTAAG".to_vec(),
                qual: Some(b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec())
            })
        );
        assert_eq!(
            rec2.unwrap().id, String::from("D00727:24:C90LYANXX:7:1109:3587:2117 1:N:0:AAGAGGCA+GTAAGGAG")
        );
        assert_eq!(
            rec3.unwrap().seq.len(), 206
        );
    }
}