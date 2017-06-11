extern crate itertools;
#[macro_use] extern crate error_chain;
#[macro_use] extern crate nom;

use std::io;

use nom::IResult;

mod internals;
use internals::*;
pub use internals::Record;
use itertools::Itertools;

error_chain! {
    errors { ParserFinished }
    foreign_links {
        Io(::std::io::Error);
    }
}

pub struct Sequences<R: io::BufRead> {
    line: Vec<u8>,
    input: R,
    parser: fn(&mut R, &mut Record, &mut Vec<u8>) -> Result<()>
}

impl<'a, R: io::BufRead> Sequences<R> {
    pub fn from_fasta(input: R) -> Sequences<R> {
        Sequences {
            line: Vec::new(),
            input, 
            parser: read_fasta2
        }
    }

    pub fn from_fastq(input: R) -> Sequences<R> {
        Sequences {
            line: Vec::new(),
            input,
            parser: read_fastq
        }
    }

    pub fn guess_format(input: R) -> Sequences<R> {
        unimplemented!()
    }

}

impl<R: io::BufRead> Iterator for Sequences<R> {
    type Item = Record;

    fn next(&mut self) -> Option<Record> {
        let mut record = Record::new();
        match (self.parser)(&mut self.input, &mut record, &mut self.line) {
            Ok(()) if record.is_empty() => None,
            Ok(()) => Some(record),
            Err(_) => None,
        }
    }
}

fn read_fasta2<R: io::BufRead>(reader: &mut R, record: &mut Record, line: &mut Vec<u8>) -> Result<()> {
    record.clear();
    if line.is_empty() {
        reader.read_until(b'\n', line)?;
        if line.is_empty() {
            return Ok(());
        }
    }

    if !line.starts_with(b">") {
        bail!("Expected > at record start.");
    }
    record.id = String::from_utf8(line[1..line.len()-1].to_vec()).unwrap();
    loop {
        line.clear();
        reader.read_until(b'\n', line)?;
        if line.is_empty() || line.starts_with(b">") {
            break;
        } else if line.ends_with(b"\n") {
            line.pop();
        }
        record.seq.append(line);
    }
    Ok(())
}

fn read_fastq<R: io::BufRead>(reader: &mut R, record: &mut Record, line: &mut Vec<u8>) -> Result<()> {
    let mut buf: Vec<u8> = Vec::new();
    loop {
        let ibuf_len: usize;
        {
            let ibuf = reader.fill_buf().expect("fastq parser: error while reading from fastq file!");  
            ibuf_len = ibuf.len();
            if ibuf_len == 0 {
                return Ok(())
            }
            buf.extend(ibuf);
        }
        match parse_single_fastq(&buf) {
            IResult::Done(extra, _record) => {
                // Consume only as much as was needed to get the record
                reader.consume(ibuf_len - extra.len());
                *record = Record{ .._record };
                return Ok(())
            },
            IResult::Error(e) => return Err("parsing error".into()),
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
        let mut seqs = Sequences::from_fasta(fa);
        let rec1 = seqs.next().unwrap();
        let rec2 = seqs.next().unwrap();
        let rec3 = seqs.next().unwrap();
        // assert_eq!(None, records.next());
        assert_eq!(
            rec1, 
            Record{
                id: String::from("gi|1348912|gb|G26680|G26680 human STS STS_D11729.^Agi|1396336|gb|G27617|G27617 human STS SHGC-32648."),
                seq: b"CGGACCAGACGGACACAGGGAGAAGCTAGTTTCTTTCATGTGATTGANATNATGACTCTACTCCTAAAAGGGAAAAANCAATATCCTTGTTTACAGAAGAGAAACAAACAAGCCCCACTCAGCTCAGTCACAGGAGAGANCACAGAAAGTCTTAGGATCATGANCTCTGAAAAAAAGAGAAACCTTATCTTTNCTTTGTGGTTCCTTTAAACACACTCACACACACTTGGTCAGAGATGCTGTGCTTCTTGGAAGCAAGGNCTCAAAGGCAAGGTGCACGCAGAGGGACGTTTGAGTCTGGGATGAAGCATGTNCGTATTATTTATATGATGGAATTTCACGTTTTTATGTNAAGCNTGACAACACCAGGCAGGTATGAGAGGAAAGCAAGGCCCGTCCATNGCTGTCCGTACNCTTACGGNTTGCTTGTNGGAGNCATTTNGGTATTGTTTGTTGTAANANCCAAAANGGGCTTTGGNNTGGNAAAAGGGCAGANNGGGGGGGTTGGTGTNGTTTTTTGGGGGGANNNTTTNGATTTGGTNCCGGGNTTTNGTTTNCCNCGGNACCGGNTTTTGGTTGGGGNCCATTTNTGNGGGGCNTTGGNGTTNCNTTNCCCNNNTNNGANTGGTTTNA"
                    .to_vec(),
                qual: None
            }
        );
        assert_eq!(
            rec2.id, String::from("gi|1348917|gb|G26685|G26685 human STS STS_D11734.")
        );
        assert_eq!(
            rec3.seq.len(), 471
        );
    }

    #[test]
    fn test_fastq_reader() {
        let fq = BufReader::new(File::open("tests/test_data/test.fq").unwrap());
        let mut records = Sequences::from_fastq(fq);
        let rec1 = records.next().unwrap();
        let rec2 = records.next().unwrap();
        let rec3 = records.next().unwrap();
        // assert_eq!(None, records.next());
        assert_eq!(
            rec1,
            Record{
                id: String::from("D00727:24:C90LYANXX:7:1109:1635:2137 1:N:0:AAGAGGCA+GTAAGGAG"),
                seq: b"GTCCAGAGCTTCGGTATAACGCTTGATCGCCAATCATTTTCGGCGCAGGATCACTCGATGAGTAAG".to_vec(),
                qual: Some(b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII".to_vec())
            }
        );
        assert_eq!(
            rec2.id, String::from("D00727:24:C90LYANXX:7:1109:3587:2117 1:N:0:AAGAGGCA+GTAAGGAG")
        );
        assert_eq!(
            rec3.seq.len(), 206
        );
    }

}