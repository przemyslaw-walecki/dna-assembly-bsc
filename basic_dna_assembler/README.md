# DNA Assembler 
## A basic DNA Assembler using de Bruijn graph

### Also contains a evaluator pipeline to compare different assemblers

## Usage:

Before using assembler/evaluator, install project:  
`poetry install`

- Assembly:  
`poetry run assemble -1 {first_read} -[2 {seconds_read}] [-k {k-mer size}] [-t {coverage threshold}] [-d {dead-end pruning}] -o {output file}`

- Evaluation:  
`poetry run compare-assemblies [--data-dir {path to assemblers outputs}] --reference {reference genome} --assemblers {space separated output names}`