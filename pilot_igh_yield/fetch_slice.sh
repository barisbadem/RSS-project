#!/usr/bin/env bash
# Download the first MiB megabytes of both mates of one ENA run (byte range, resumable).
#
# usage: fetch_slice.sh RUN OUTDIR [MiB=300]
#
# FASTQ reads are in flowcell order, so the head of the file is a random sample of
# the library. A 300 MiB slice of a 2x75 bp run is roughly 9 million read pairs.
# The ENA tunnel can drop a long transfer; the loop resumes from the bytes already
# on disk, so re-running the script after a failure is safe.
set -euo pipefail

run="$1"; out="$2"; mib="${3:-300}"
mkdir -p "$out/$run"
want=$(( mib * 1024 * 1024 ))

# filereport prepends the run accession as column 1; the paired URLs are in column 2, ';'-separated.
urls=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=fastq_ftp&format=tsv&limit=0" | tail -n 1 | cut -f2)
u1=$(tr ';' '\n' <<< "$urls" | grep '_1\.fastq\.gz$' | head -n 1 || true)
u2=$(tr ';' '\n' <<< "$urls" | grep '_2\.fastq\.gz$' | head -n 1 || true)
if [ -z "$u1" ] || [ -z "$u2" ]; then
  echo "no paired FASTQ URLs found for ${run}: ${urls}" >&2
  exit 1
fi

for mate in 1 2; do
  u=$([ "$mate" = 1 ] && echo "$u1" || echo "$u2")
  f="$out/$run/R${mate}.partial.gz"
  cur=$(stat -c %s "$f" 2>/dev/null || echo 0)
  tries=0
  while [ "$cur" -lt "$want" ] && [ "$tries" -lt 8 ]; do
    curl -sS --max-time 1500 -r "${cur}-$(( want - 1 ))" "https://${u}" >> "$f" || sleep 3
    cur=$(stat -c %s "$f")
    tries=$(( tries + 1 ))
  done
  echo "$run R${mate}: ${cur} bytes"
done
