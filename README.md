# Fastx

A library for reading FASTA and FASTQ files in Rust. 

## Advantages

- Returns generic `Record` structs that are input-agnostic, i.e. not typed by the input file type, making it suitable for applications that need to handle either FASTA or FASTQ files.
- Extremely fast parsing using [nom](http://rust.unhandledexpression.com/nom/index.html); can handle input from stdin or from file equally well.

## Usage

```rust
extern crate fastx;
use std::io::{File, BufReader};
use fastx::Records;

let input = BufReader::new(File::open("foo.txt"));
let records = Records::from_fasta(input);
for record in records {
    println!("{}", record.id);
}
```

