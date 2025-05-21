#!/bin/bash

poetry run compare-assemblies --data-dir ./data --reference ./data/GCF_000007365.1_ASM736v1_genomic.fna --assemblers dnaasm --output report.csv --md-log comparison.log.md
