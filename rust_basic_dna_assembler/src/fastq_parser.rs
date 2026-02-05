//! Moduł parsera FASTQ.
//!
//! Udostępnia minimalny parser służący do odczytu sekwencji DNA z plików
//! FASTQ w standardowym formacie czteroliniowym.  
//! Parser pomija identyfikatory i jakość — zwraca wyłącznie surowe sekwencje DNA.

use std::fs::File;
use std::io::{self, BufRead, BufReader};

/// Prosty parser plików FASTQ.
///
/// Odczytuje jedynie linie sekwencji (druga linia każdego rekordu FASTQ).
pub struct FastqParser {
    /// Ścieżka do pliku FASTQ.
    pub path: String,
}

impl FastqParser {
    /// Tworzy nowy parser na podstawie ścieżki do pliku.
    ///
    /// # Argumenty
    /// * `path` — ścieżka do pliku FASTQ.
    pub fn new(path: impl Into<String>) -> Self {
        Self { path: path.into() }
    }

    /// Wczytuje sekwencje DNA z pliku FASTQ.
    ///
    /// Zwraca:
    /// - `Vec<String>`: lista sekwencji DNA.
    ///
    /// # Błędy
    /// Zwraca `io::Error`, jeśli plik nie może zostać otwarty lub odczytany.
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

            rdr.read_line(&mut line)?; // linia sekwencji
            reads.push(line.trim_end().to_string());
            line.clear();

            rdr.read_line(&mut line)?; // linia "+"
            line.clear();

            rdr.read_line(&mut line)?; // linia jakości
            line.clear();
        }

        Ok(reads)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn tmp_fastq_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("fastq_parser_test_{}_{}_{}.fastq", name, std::process::id(), ts));
        p
    }

    fn write_text(path: &PathBuf, content: &str) {
        let mut f = fs::File::create(path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f.flush().unwrap();
    }

    #[test]
    fn load_reads_parses_two_records() {
        let fastq = "\
@r1
ACGT
+
!!!!
@r2
TTAA
+
####\n";

        let path = tmp_fastq_path("two_records");
        write_text(&path, fastq);

        let parser = FastqParser::new(path.to_string_lossy().to_string());
        let reads = parser.load_reads().unwrap();

        fs::remove_file(&path).unwrap();

        assert_eq!(reads, vec!["ACGT".to_string(), "TTAA".to_string()]);
    }

    #[test]
    fn load_reads_trims_newline() {
        let fastq = "\
@r1
ACGT
+
!!!!\n";
    
        let path = tmp_fastq_path("trim_newline");
        write_text(&path, fastq);
    
        let parser = FastqParser::new(path.to_string_lossy().to_string());
        let reads = parser.load_reads().unwrap();
    
        fs::remove_file(&path).unwrap();
    
        assert_eq!(reads, vec!["ACGT".to_string()]);
        assert_eq!(reads[0].as_bytes().last().copied(), Some(b'T'));
    }

    #[test]
    fn load_reads_empty_file_returns_empty_vec() {
        let path = tmp_fastq_path("empty");
        write_text(&path, "");

        let parser = FastqParser::new(path.to_string_lossy().to_string());
        let reads = parser.load_reads().unwrap();

        fs::remove_file(&path).unwrap();

        assert!(reads.is_empty());
    }

    #[test]
    fn load_reads_nonexistent_file_returns_error() {
        let path = tmp_fastq_path("missing");
        // nie tworzymy pliku

        let parser = FastqParser::new(path.to_string_lossy().to_string());
        let err = parser.load_reads().unwrap_err();

        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
    }
}
