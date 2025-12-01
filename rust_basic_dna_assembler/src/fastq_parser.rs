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
