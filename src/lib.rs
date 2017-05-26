#[macro_use]
extern crate nom;
extern crate regex;
extern crate itertools;

use std::io;
use std::str;

use itertools::Itertools;
use nom::{IResult, not_line_ending, line_ending, multispace, space, rest};

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Record<'a> {
    pub id: &'a str,
    pub seq: Vec<&'a u8>,
    pub qual: Option<&'a [u8]>
}

named!(parse_fasta<&[u8], Vec<Record>>,
    many1!(do_parse!(
        tag!(">") >>
        id: map_res!(not_line_ending, str::from_utf8) >> line_ending >>
        seq: is_not!(">") >>
        (Record{ 
            id: id.trim(), 
            seq: seq.iter().filter(|b| **b != '\n' as u8).collect::<Vec<&u8>>(),
            qual: None,
        })
    ))
);

named!(quality_scores, 
    do_parse!(
        tag!("+") >> line_ending >>
        scores: not_line_ending >>
        (scores)
    )
);

named!(parse_fastq<&[u8], Vec<Record>>,
    many1!(dbg_dmp!(do_parse!(
        tag!("@") >>
        id: map_res!(not_line_ending, str::from_utf8) >> line_ending >>
        seq: not_line_ending >> line_ending >>
        qual: opt!(quality_scores) >> line_ending >>
        (Record{ id: id.trim(), qual, seq: seq.iter().filter(|b| **b != '\n' as u8).collect::<Vec<&u8>>() })
    )))
);

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_fasta() {
        let bytes = b">id1 desc \nTTTTTT\nAAAAAA\n>id2\nGGGGGG\nCCCCCC\n";
        let records = parse_fasta(bytes).unwrap().1;
        let rec1 = Record{id:"id1 desc", seq: b"TTTTTTAAAAAA".iter().collect(), qual: None};
        let rec2 = Record{id:"id2", seq: b"GGGGGGCCCCCC".iter().collect(), qual: None};
        assert_eq!(rec1, records[0]);
        assert_eq!(rec2, records[1]);
    }

    #[test]
    fn test_parse_fastq() {
        let bytes = b"@id1 desc \nTTTTTTAAAAAA\n+\nFFFFFFFFFFFFF\n@id2\nGGGGGGCCCCCC\n+\nFFFFFFFFFFFFF\n";
        let records = parse_fastq(bytes).unwrap().1;
        let rec1 = Record{id:"id1 desc", seq: b"TTTTTTAAAAAA".iter().collect(), qual: Some(b"FFFFFFFFFFFFF")};
        let rec2 = Record{id:"id2", seq: b"GGGGGGCCCCCC".iter().collect(), qual: Some(b"FFFFFFFFFFFFF")};
        assert_eq!(rec1, records[0]);
        assert_eq!(rec2, records[1]);
    }
}
