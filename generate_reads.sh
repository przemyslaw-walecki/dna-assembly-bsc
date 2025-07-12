docker run --rm \
  -v "$(pwd)/random_genome_generator:/data" \
  -w /data \
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 \
  sh -c "pirs simulate random_50.fna -x 200 -l 35 -m 35 -v 2 -e 0.00 -o reads && mv reads_1.fq reads50_v2.fq"
