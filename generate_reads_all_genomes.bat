@echo off
setlocal enabledelayedexpansion

REM ---- Directory with genomes (.fasta) ----
set GENOME_DIR="C:\Users\Przemek\Desktop\Inzynierka stuff\own_assembler\dna-assembly-bsc\genomes"

REM ---- Change to genome directory ----
cd %GENOME_DIR%

REM ---- For each FASTA file ----
for %%F in (*.fasta) do (

    REM Full filename: e.g. Klebsiella-pneu.fasta
    set FILE=%%F

    REM Extract first token before "-" 
    for /f "tokens=1 delims=-" %%A in ("%%F") do (
        set PREFIX=%%A
    )

    echo ===================================================
    echo Processing genome: %%F
    echo Prefix: !PREFIX!
    echo ===================================================

    REM Run pirs in Docker
    docker run --rm ^
        -v %GENOME_DIR%:/data ^
        -w /data ^
        biocontainers/pirs:v2.0.2dfsg-8-deb_cv1 ^
        sh -c "pirs simulate %%F -x 20 -l 100 -m 300 -v 2 -e 0.00 -o !PREFIX!_clean_cov20"

    echo.
)

echo DONE!
pause
