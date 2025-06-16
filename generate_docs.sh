#!/bin/bash
set -e

# Generowanie dokumentacji Rust
cd rust_basic_dna_assembler
cargo doc --no-deps
cd ..

# Kopiowanie wygenerowanej dokumentacji Rust do katalogu MkDocs
mkdir -p docs/static/rustdoc
cp -r rust_basic_dna_assembler/target/doc/* docs/static/rustdoc/

# Generowanie dokumentacji MkDocs (Python + statyczne strony)
mkdocs build
