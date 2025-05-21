use std::fs::File;
use std::io::{self, BufRead, BufReader};

pub struct FastqParser {
    pub path: String,
}

impl FastqParser {
    pub fn new(path: impl Into<String>) -> Self {
        Self { path: path.into() }
    }

    pub fn load_reads(&self) -> io::Result<Vec<String>> {
        let file = File::open(&self.path)?;
        let mut rdr = BufReader::new(file);
        let mut reads = Vec::new();
        let mut line = String::new();

        loop {

            if rdr.read_line(&mut line)? == 0 {
                break;
            }
            line.clear();

            rdr.read_line(&mut line)?;
            reads.push(line.trim_end().to_string());
            line.clear();

            rdr.read_line(&mut line)?;
            line.clear();

            rdr.read_line(&mut line)?;
            line.clear();
        }

        Ok(reads)
    }
}
