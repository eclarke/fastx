use std::str;
use itertools::Itertools;

use nom::{not_line_ending, line_ending, rest};

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct Record {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>
}

named!(pub parse_fasta_header<&[u8], &str>,
    do_parse!(
        tag!(">") >>
        id: map_res!(not_line_ending, str::from_utf8) >>
        line_ending >>
        (id))
);

named!(pub parse_fasta_body<&[u8], Vec<char>>,
    do_parse!(
        body: many1!(ws!(none_of!(">"))) >>
        peek!(tag!(">")) >>
        (body))
);

named!(pub parse_single_fasta<&[u8], Record>,
    do_parse!(
        id: parse_fasta_header >>
        seq: parse_fasta_body >>
        (Record {
            id: String::from(id.trim()),
            seq: seq.into_iter().map(|c| c as u8).collect::<Vec<u8>>(),
            qual: None
        })
    )
);

named!(pub parse_single_fasta1<&[u8], Record>,
    do_parse!(
        tag!(">") >>
        id: map_res!(not_line_ending, str::from_utf8) >>
        // seq: map!(take_until!(">"), deref) >>
        seq: many_till!(
                do_parse!(
                    s: map!(not_line_ending, deref) >>
                    line_ending >>
                    (s)), peek!(tag!(">"))) >>
        (Record{ 
            id: String::from(id.trim()),
            seq: seq.0.into_iter().flatten().collect::<Vec<u8>>(),
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
    fn test_parse_fa_header() {
        let bytes = b">id\nBLAHBLABAH";
        let (rest, id) = parse_fasta_header(bytes).unwrap();
        assert_eq!("id", id);
    }

    #[test]
    fn test_parse_fa_body() {
        let bytes = b"AA\nTT\nGG\n>skip";
        let (rest, body) = parse_fasta_body(bytes).unwrap();
        assert_eq!("AATTGG", body.into_iter().collect::<String>());
        assert_eq!(b">skip".to_vec(), rest);
    }

    // #[test]
    // fn test_parse_fa_body_noend() {
    //     let bytes = b"AA\nTT\nGG\n";
    //     let (rest, body) = parse_fasta_body(bytes).unwrap();
    //     assert_eq!("AATTGG", body.into_iter().collect::<String>());
    // }

    // #[test]
    // fn test_parse_fasta() {
    //     let bytes = b">id1 desc \nTTTTTT\nAAAAAA\n>id2 desc\nBLAH\n";
    //     let (rest, record1) = parse_single_fasta(bytes).unwrap();
    //     let rec1 = Record{id:String::from("id1 desc"), seq: b"TTTTTTAAAAAA".to_vec(), qual: None};
    //     assert_eq!(rec1, record1);
    //     assert_eq!(b">id2 desc\nBLAH\n", rest);
    // }

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