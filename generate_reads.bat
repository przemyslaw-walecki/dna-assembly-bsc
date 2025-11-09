docker run --rm ^
  -v "C:\Users\Przemek\Desktop\Inzynierka stuff\own_assembler\dna-assembly-bsc\genomes\:/data" ^
  -w /data ^
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 ^
  sh -c "pirs simulate Escherichia_coli_str_K-12_substr_MG1655.fasta -x 10 -l 100 -m 300 -v 2 -e 0.00 -o ecoli_clean_cov10"
