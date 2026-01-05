#!/usr/bin/env bash
set -euo pipefail

# =========================
# Config (repo-relative)
# =========================
K_FIXED=21

ASSEMBLER_BIN="rust_basic_dna_assembler/target/release/rust_basic_dna_assembler"
GNN_SCRIPT="bubble_resolution_gnn/GNN_infer.py"

# Default locations (recommended to create/move)
DEFAULT_MODEL_DIR="models/default"
DEFAULT_CKPT="${DEFAULT_MODEL_DIR}/checkpoint.pth"

# Defaults for external tools/interpreters
PYTHON_BIN="${PYTHON_BIN:-python3}"
BUBBLEGUN_BIN="${BUBBLEGUN_BIN:-BubbleGun}"

# =========================
# Helpers
# =========================
usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh --reads-1 <R1.fq> --reads-2 <R2.fq> --out <out_dir> [--sample <name>] [--mode gnn|coverage|none]
                   [--model-dir <dir> | --ckpt <path>] [--cpu]

Notes:
- k-mer size is FIXED to k=21 and cannot be changed (models trained only for k=21).
- Default mode is: gnn
- Default model checkpoint: models/default/checkpoint.pth

Examples:
  ./run_pipeline.sh --reads-1 genomes/sample_1.fq --reads-2 genomes/sample_2.fq --out out --sample salmonella --mode gnn --cpu
  ./run_pipeline.sh --reads-1 genomes/sample_1.fq --reads-2 genomes/sample_2.fq --out out --sample salmonella --mode coverage
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }

need_file() { [[ -f "$1" ]] || die "Missing file: $1"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Missing command in PATH: $1"; }

# =========================
# Args
# =========================
READS_1=""
READS_2=""
OUT_ROOT=""
SAMPLE=""
MODE="gnn"
MODEL_DIR=""
CKPT=""
USE_CPU="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads-1|-1) READS_1="${2:-}"; shift 2;;
    --reads-2|-2) READS_2="${2:-}"; shift 2;;
    --out|-o) OUT_ROOT="${2:-}"; shift 2;;
    --sample) SAMPLE="${2:-}"; shift 2;;
    --mode) MODE="${2:-}"; shift 2;;
    --model-dir) MODEL_DIR="${2:-}"; shift 2;;
    --ckpt) CKPT="${2:-}"; shift 2;;
    --cpu) USE_CPU="true"; shift 1;;
    --k|-k)
      # Hard refuse to prevent silent mismatch
      die "Parameter k is not configurable. This pipeline supports only k=21."
      ;;
    -h|--help) usage; exit 0;;
    *) die "Unknown argument: $1 (use --help)";;
  esac
done

[[ -n "$READS_1" && -n "$READS_2" && -n "$OUT_ROOT" ]] || { usage; exit 1; }

case "$MODE" in
  gnn|coverage|none) ;;
  *) die "Invalid --mode: $MODE (expected: gnn|coverage|none)";;
esac

# =========================
# Validate environment
# =========================
need_file "$READS_1"
need_file "$READS_2"
need_cmd "$PYTHON_BIN"

# BubbleGun needed if we are not skipping it inside assembler; your assembler calls it via --bubblegun-bin
# If BubbleGun is not in PATH, user can export BUBBLEGUN_BIN=/path/to/BubbleGun
need_cmd "$BUBBLEGUN_BIN"

[[ -x "$ASSEMBLER_BIN" ]] || die "Assembler binary not found/executable: $ASSEMBLER_BIN
Build it first:
  (cd rust_basic_dna_assembler && cargo build --release)
"

[[ -f "$GNN_SCRIPT" ]] || die "GNN script not found: $GNN_SCRIPT"

# Resolve checkpoint
if [[ -n "$MODEL_DIR" && -n "$CKPT" ]]; then
  die "Use either --model-dir or --ckpt, not both."
fi

if [[ -n "$MODEL_DIR" ]]; then
  # Convention: checkpoint file inside model dir
  CKPT="${MODEL_DIR%/}/checkpoint.pth"
fi

if [[ -z "$CKPT" ]]; then
  CKPT="$DEFAULT_CKPT"
fi

if [[ "$MODE" == "gnn" ]]; then
  need_file "$CKPT"
fi

# sample name: if not provided, infer from reads_1 filename (strip extension)
if [[ -z "$SAMPLE" ]]; then
  base="$(basename "$READS_1")"
  SAMPLE="${base%%.*}"
fi

# Output layout
OUT_DIR="${OUT_ROOT%/}/${SAMPLE}"
WORK_DIR="${OUT_DIR}/work"
mkdir -p "$WORK_DIR"

GFA_OUT="${WORK_DIR}/${SAMPLE}_k${K_FIXED}_graph.gfa"
BUBBLES_JSONL="${WORK_DIR}/${SAMPLE}_k${K_FIXED}_bubbles.jsonl"
DECISIONS_JSONL="${WORK_DIR}/${SAMPLE}_k${K_FIXED}_decisions.jsonl"
CONTIGS_OUT="${OUT_DIR}/${SAMPLE}_${MODE}_k${K_FIXED}.fa"

# =========================
# Build GNN args (passed as a single string to assembler)
# =========================
GNN_ARGS="--gfa ${GFA_OUT} --ckpt ${CKPT}"
if [[ "$USE_CPU" == "true" ]]; then
  GNN_ARGS="${GNN_ARGS} --cpu"
fi

# =========================
# Run
# =========================
echo "[run_pipeline] mode=${MODE}, k=${K_FIXED}"
echo "[run_pipeline] reads: ${READS_1} ${READS_2}"
echo "[run_pipeline] out:   ${CONTIGS_OUT}"
echo "[run_pipeline] work:  ${WORK_DIR}"

case "$MODE" in
  gnn)
    "$ASSEMBLER_BIN" \
      -1 "$READS_1" \
      -2 "$READS_2" \
      -k "$K_FIXED" \
      -o "$CONTIGS_OUT" \
      --gfa-out "$GFA_OUT" \
      --bubblegun-bin "$BUBBLEGUN_BIN" \
      --bubbles-jsonl "$BUBBLES_JSONL" \
      --python-bin "$PYTHON_BIN" \
      --gnn-script "$GNN_SCRIPT" \
      --decisions-jsonl "$DECISIONS_JSONL" \
      --gnn-args="$GNN_ARGS"
    ;;
  coverage)
    "$ASSEMBLER_BIN" \
      -1 "$READS_1" \
      -2 "$READS_2" \
      -k "$K_FIXED" \
      -o "$CONTIGS_OUT" \
      --bubble-resolver coverage \
      --bubblegun-bin "$BUBBLEGUN_BIN" \
      --gfa-out "$GFA_OUT" \
      --bubbles-jsonl "$BUBBLES_JSONL" \
      --skip-gnn
    ;;
  none)
    "$ASSEMBLER_BIN" \
      -1 "$READS_1" \
      -2 "$READS_2" \
      -k "$K_FIXED" \
      -o "$CONTIGS_OUT" \
      --bubble-resolver none \
      --bubblegun-bin "$BUBBLEGUN_BIN" \
      --gfa-out "$GFA_OUT" \
      --bubbles-jsonl "$BUBBLES_JSONL" \
      --skip-gnn
    ;;
esac

echo "[run_pipeline] done"
echo "[run_pipeline] contigs:   $CONTIGS_OUT"
echo "[run_pipeline] graph:     $GFA_OUT"
echo "[run_pipeline] bubbles:   $BUBBLES_JSONL"
if [[ "$MODE" == "gnn" ]]; then
  echo "[run_pipeline] decisions: $DECISIONS_JSONL"
fi
