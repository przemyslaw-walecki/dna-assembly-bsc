#!/bin/bash
set -e

echo "[1] Generowanie dokumentacji Rust..."

cd rust_basic_dna_assembler
cargo doc --no-deps
cd ..

cd rust_bubble_labelling
cargo doc --no-deps
cd ..

echo "[2] Kopiowanie dokumentacji Rust..."

# Usuwamy stare rzeczy
rm -rf docs/static/rustdoc
mkdir -p docs/static/rustdoc/basic
mkdir -p docs/static/rustdoc/bubble

# Kopiujemy assembler
cp -r rust_basic_dna_assembler/target/doc/* docs/static/rustdoc/basic/

# Kopiujemy bubble-labelling
cp -r rust_bubble_labelling/target/doc/* docs/static/rustdoc/bubble/

# Czyścimy licencje Adobe
find docs/static/rustdoc -name "SourceSerif4-LICENSE-*" -delete || true

echo "[3] Budowanie MkDocs..."
mkdocs build

echo "[DONE]"
