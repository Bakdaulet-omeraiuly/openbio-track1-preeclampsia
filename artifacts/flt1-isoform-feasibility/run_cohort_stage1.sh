#!/bin/bash
# Stage 1 of Addendum J1: fetch one SRA run, dump, filter, verify read count, delete the download.
# Usage: run_cohort_stage1.sh WORKDIR RUN EXPECTED_SPOTS     (resumable: skips runs already verified)
set -u
WORK="$1"; RUN="$2"; SPOTS="$3"
cd "$WORK" || exit 2
[ -s "logs/$RUN.ok" ] && exit 0
for attempt in 1 2 3; do
  # keep at least ~2 GB free before starting another download
  while [ "$(df -k . | awk 'NR==2{print $4}')" -lt 2000000 ]; do sleep 20; done
  rm -f "sra/$RUN.sra" "cand/$RUN.fastq" "cand/$RUN.fastq.gz"
  if curl -sf --retry 3 -o "sra/$RUN.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN" &&
     env/bin/fastq-dump -Z "sra/$RUN.sra" 2>"logs/$RUN.dump.err" |
       ./flt1_stream_filter ref/marker-kmers.tsv ref/flt1-bait.fa "$RUN" "cand/$RUN.fastq" "logs/$RUN.json"; then
    got=$(grep -o '"total_reads": [0-9]*' "logs/$RUN.json" | grep -o '[0-9]*$')
    bad=$(grep -o '"invalid_fastq_records": [0-9]*' "logs/$RUN.json" | grep -o '[0-9]*$')
    if [ "$got" = "$SPOTS" ] && [ "$bad" = "0" ]; then
      gzip -f "cand/$RUN.fastq" && rm -f "sra/$RUN.sra"
      echo "$RUN $got $(date -u +%FT%TZ) attempt=$attempt" > "logs/$RUN.ok"
      exit 0
    fi
    echo "$RUN read-count mismatch: got=$got expected=$SPOTS invalid=$bad (attempt $attempt)" >> logs/stage1.errors
  else
    echo "$RUN fetch/dump failed (attempt $attempt)" >> logs/stage1.errors
  fi
  rm -f "sra/$RUN.sra"
done
exit 1
