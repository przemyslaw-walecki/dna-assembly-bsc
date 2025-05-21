@echo off
if not exist "%~dp0config.bat" (
  echo ERROR: config.bat not found! Copy config.bat.example and fill it in.
  exit /b 1
)
call "%~dp0config.bat"

if "%DATA_DIR%"=="" (
  echo ERROR: DATA_DIR not set in config.bat
  exit /b 1
)
if "%GENOME_LENGTH%"=="" (
  echo ERROR: GENOME_LENGTH not set in config.bat
  exit /b 1
)
docker run --rm ^
  -v "%DATA_DIR%:/data" ^
  -w /data ^
  wkusmirek/dnaasm:latest ^
  dnaasm -assembly ^
    -k 55 ^
    -genome_length "%GENOME_LENGTH%" ^
    -paired_reads_algorithm 1 ^
    -insert_size_mean_inward 400 ^
    -insert_size_std_dev_inward 20 ^
    -single_edge_counter_threshold 5 ^
    -i1_1 "%READ_1%" ^
    -i1_2 "%READ_2%" ^
    -output_file_name "%OUTPUT_FILE%"