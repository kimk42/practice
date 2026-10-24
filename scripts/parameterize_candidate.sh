#!/usr/bin/env bash
# GAFF2 + AM1-BCC exploratory assignment; retain full logs and fail on tool errors.
set -euo pipefail
project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$project_dir/scripts/activate.sh"
if [ "$#" -ne 3 ]; then
    echo 'Usage: parameterize_candidate.sh INPUT_MOL NAME NEW_OUTPUT_DIR' >&2
    exit 2
fi
input=$(realpath -- "$1")
name=$2
output=$(realpath -m -- "$3")
if [ -e "$output/$name.acpype" ]; then
    echo "Output already exists: $output/$name.acpype; preserve it and use a fresh directory" >&2
    exit 2
fi
mkdir -p "$output"
export LD_LIBRARY_PATH="$MD_TOOLS_DIR/root/usr/lib/x86_64-linux-gnu/openblas-pthread${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# Discover bundled tools explicitly so ACPYPE does not reset LD_LIBRARY_PATH.
export AMBERHOME=$(python -c 'from acpype.utils import bundled_amber_dir; print(bundled_amber_dir())')
export PATH="$AMBERHOME/bin:$PATH"
ldd "$AMBERHOME/bin/sqm" > "$output/runtime-libraries.log"
if ! rg -q 'openblas-pthread' "$output/runtime-libraries.log"; then
    echo 'Optimized BLAS was not resolved; stopping before the expensive calculation' >&2
    exit 2
fi
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd "$output"
acpype -i "$input" -b "$name" -a gaff2 -c bcc -n 0 -m 1 -s "${MD_SQM_MAX_SECONDS:-10800}" > "$name.log" 2>&1
