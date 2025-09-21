docker run --rm ^
  -v "C:\Users\Przemek\Desktop\Inzynierka stuff\own_assembler\dna-assembly-bsc\genomes\:/data" ^
  -w /data ^
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 ^
  sh -c "pirs simulate Klebsiella_pneumoniae_subsp_pneumoniae_NTUH-K2044_DNA.fasta -x 20 -l 100 -m 300 -v 2 -e 0.00 -o klebsiella_clean_cov20"
