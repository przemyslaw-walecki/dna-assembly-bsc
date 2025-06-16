docker run --rm `
  -v "$($(pwd).Path)\random_genome_generator:/data" `
  -w /data `
  biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 `
  pirs simulate random.fna -x 50 -l 100 -m 400 -v 20 -e 0.01 -o reads
