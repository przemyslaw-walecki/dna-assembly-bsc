//! FASTQ parser module.
//!
//! Provides a minimal parser for reading sequencing reads from a FASTQ file.
//! It supports extracting DNA sequences from files encoded in the standard
//! four-line FASTQ format.
//!
//! # Overview
//!
//! The `FastqParser` struct exposes a `load_reads` method which reads the file,
//! skips quality and identifier lines, and collects only the sequence lines.
//!
//! # Limitations
//!
//! - Assumes the FASTQ format is strictly correct (exactly four lines per record).
//! - Reads are returned as strings without additional metadata or quality scores.

use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// A simple parser for FASTQ files.
///
/// The parser reads only the sequence lines from the standard FASTQ format.
/// It does not parse quality scores or read identifiers.
pub struct FastqParser {
    /// Path to the FASTQ file.
    pub path: String,
}

impl FastqParser {
    /// Creates a new `FastqParser` from a given file path.
    ///
    /// # Arguments
    /// * `path` - Path to the FASTQ file (as `String` or `&str`)
    pub fn new(path: impl Into<String>) -> Self {
        Self { path: path.into() }
    }

    
    /// Loads and parses DNA reads from the FASTQ file.
    ///
    /// Only the second line of each FASTQ record (the sequence) is extracted.
    ///
    /// # Returns
    /// A vector of strings, each representing a DNA sequence read.
    ///
    /// # Errors
    /// Returns an `io::Error` if the file cannot be opened or read.
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

            rdr.read_line(&mut line)?;  // sequence line
            reads.push(line.trim_end().to_string());
            line.clear();

            rdr.read_line(&mut line)?;  // plus line
            line.clear();

            rdr.read_line(&mut line)?;  // quality line
            line.clear();
        }

        Ok(reads)
    }
}
