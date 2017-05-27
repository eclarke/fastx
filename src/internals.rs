use std::str;
use itertools::Itertools;

use nom::{not_line_ending, line_ending};

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>
}

named!(pub parse_single_fasta<&[u8], Record>,
    do_parse!(
        tag!(">") >>
        id: map_res!(not_line_ending, str::from_utf8) >> line_ending >>
        seq: many1!(do_parse!(
            s: map!(not_line_ending, deref) >> 
            line_ending >> (s)
        )) >> 
        (Record{ 
            id: String::from(id.trim()), 
            seq: seq.into_iter().flatten().collect::<Vec<u8>>(),
            qual: None,
        })
    )
);

fn deref(s: &[u8]) -> Vec<u8> {
    s.into_iter().map(|b| *b).collect::<Vec<u8>>()
}

named!(quality_scores<&[u8], Vec<u8>>, 
    do_parse!(
        tag!("+") >> line_ending >>
        scores: map!(not_line_ending, deref) >>
        (scores)
    )
);

named!(pub parse_single_fastq<&[u8], Record>,
    do_parse!(
        tag!("@") >>
        id: map_res!(not_line_ending, str::from_utf8) >> line_ending >>
        seq: map!(not_line_ending, deref) >> line_ending >>
        qual: opt!(quality_scores) >> line_ending >>
        (Record{ 
            id: String::from(id.trim()), 
            qual: qual,
            seq: seq })
    )
);

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_parse_fasta() {
        let bytes = b">id1 desc \nTTTTTT\nAAAAAA\n>id2\nGGGGGG\nCCCCCC\n";
        let (rest, record1) = parse_single_fasta(bytes).unwrap();
        let (_, record2) = parse_single_fasta(rest).unwrap();
        let rec1 = Record{id:String::from("id1 desc"), seq: b"TTTTTTAAAAAA".to_vec(), qual: None};
        let rec2 = Record{id:String::from("id2"), seq: b"GGGGGGCCCCCC".to_vec(), qual: None};
        assert_eq!(rec1, record1);
        assert_eq!(rec2, record2);
    }

    #[test]
    fn test_parse_fastq() {
        let bytes = b"@id1 desc \nTTTTTTAAAAAA\n+\nFFFFFFFFFFFFF\n@id2\nGGGGGGCCCCCC\n+\nFFFFFFFFFFFFF\n";
        let (rest, record1) = parse_single_fastq(bytes).unwrap();
        let (_, record2) = parse_single_fastq(rest).unwrap();
        let rec1 = Record{id:String::from("id1 desc"), seq: b"TTTTTTAAAAAA".to_vec(), qual: Some(b"FFFFFFFFFFFFF".to_vec())};
        let rec2 = Record{id:String::from("id2"), seq: b"GGGGGGCCCCCC".to_vec(), qual: Some(b"FFFFFFFFFFFFF".to_vec())};
        assert_eq!(rec1, record1);
        assert_eq!(rec2, record2);
    }
}