docker run --rm ^
  -v "C:\Users\Przemek\Desktop\Inzynierka stuff\own_assembler\dna-assembly-bsc\genomes\:/data" ^
  -w /data ^
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 ^
  sh -c "pirs simulate Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_LT2.fasta -x 50 -l 100 -m 300 -v 2 -e 0.00 -o salmonella_clean && mv salmonella_clean_1.fq salmonella_clean_2.fq"
