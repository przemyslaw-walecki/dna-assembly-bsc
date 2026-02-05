docker run --rm ^
  -v "C:\Users\Przemek\Desktop\Inzynierka stuff\own_assembler\dna-assembly-bsc\genomes\ecoli_variants\:/data" ^
  -w /data ^
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 ^
  sh -c "pirs simulate Ecoli_variant5.fasta -x 20 -l 100 -m 300 -v 2 -e 0.00 -o ecoli_variant5_clean_cov20"
