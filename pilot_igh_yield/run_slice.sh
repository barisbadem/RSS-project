#!/usr/bin/env bash
# Run TRUST4 (mouse IMGT reference) on one downloaded slice, then delete the plain FASTQ.
#
# usage: TRUST4_DIR=/path/to/TRUST4 run_slice.sh SLICE_DIR OUT_DIR RUN
#
# SLICE_DIR holds R1.partial.gz and R2.partial.gz (see fetch_slice.sh).
# TRUST4_DIR is a checkout of github.com/liulab-dfci/TRUST4 on which `make` has run;
# its mouse/ folder carries mouse_IMGT+C.fa. Env: THREADS (default 4).
set -euo pipefail

slice="$1"; out="$2"; run="$3"
t4="${TRUST4_DIR:?set TRUST4_DIR to the TRUST4 checkout (make already run)}"
here="$(cd "$(dirname "$0")" && pwd)"

mkdir -p "$out/$run"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

pairs=$(python3 -I "$here/slice_pair.py" "$slice/R1.partial.gz" "$slice/R2.partial.gz" "$tmp/R1.fq" "$tmp/R2.fq")
echo "pairs_used=${pairs}" > "$out/$run/pairs_used.txt"

"$t4/run-trust4" -f "$t4/mouse/mouse_IMGT+C.fa" --ref "$t4/mouse/mouse_IMGT+C.fa" \
    -1 "$tmp/R1.fq" -2 "$tmp/R2.fq" -t "${THREADS:-4}" \
    -o "$run" --od "$out/$run" --outputReadAssignment > "$out/$run/trust4.log" 2>&1

echo "trust4_done=1" >> "$out/$run/pairs_used.txt"
