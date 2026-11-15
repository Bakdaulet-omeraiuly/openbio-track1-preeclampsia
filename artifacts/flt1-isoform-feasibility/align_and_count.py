"""Stage 2 of Addendum J1: HISAT2 full-genome alignment of FLT1 candidate reads and marker counting.

Counting rules are those fixed in raw-reanalysis-plan.md, Addendum J1, item 4. No clinical group label
is read by this script.

Usage:
  python align_and_count.py --work WORK --markers            # synthetic marker 50-mer confirmation
  python align_and_count.py --work WORK --runs SRR1 SRR2 ...  # donors (default: every cand/*.fastq.gz)
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
FLT1_SPAN = ("13", 28874489, 29069265)
PRIMARY_FLANK = 10
SENSITIVITY_FLANK = 5


def load_markers():
    with (HERE / "flt1-marker-definitions.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    junctions = {}
    bodies = {}
    for row in rows:
        key = (int(row["start_grch37"]), int(row["end_grch37"]))
        if row["marker_type"] == "unique_splice_junction_within_flt1":
            junctions[key] = row["marker_id"]
        else:
            bodies[row["marker_id"]] = key
    return [row["marker_id"] for row in rows], junctions, bodies


def parse_alignment(pos, cigar):
    """Return (reference end, introns [(start, end, aligned_left, aligned_right)], spliced?)."""
    ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    aligned_total = sum(n for n, op in ops if op in "M=X")
    ref = pos
    aligned_before = 0
    introns = []
    for n, op in ops:
        if op in "M=X":
            ref += n
            aligned_before += n
        elif op == "D":
            ref += n
        elif op == "N":
            introns.append((ref, ref + n - 1, aligned_before, aligned_total - aligned_before))
            ref += n
    return ref - 1, introns


def count_sam(lines, marker_ids, junctions, bodies, name_hook=None):
    counts = {m: 0 for m in marker_ids}
    flank5 = {m: 0 for m in marker_ids if "junction" in m}
    stats = {"aligned_records": 0, "unique_reads": 0, "unique_reads_in_flt1_span": 0,
             "flt1_span_forward_strand": 0, "flt1_span_reverse_strand": 0}
    per_read = {}
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        flag = int(fields[1])
        if flag & 4 or flag & 256 or flag & 2048:
            continue
        stats["aligned_records"] += 1
        nh = next((int(tag[5:]) for tag in fields[11:] if tag.startswith("NH:i:")), None)
        if nh != 1:
            continue
        stats["unique_reads"] += 1
        chrom, pos, cigar = fields[2], int(fields[3]), fields[5]
        if chrom != FLT1_SPAN[0]:
            continue
        end, introns = parse_alignment(pos, cigar)
        if end < FLT1_SPAN[1] or pos > FLT1_SPAN[2]:
            continue
        stats["unique_reads_in_flt1_span"] += 1
        stats["flt1_span_reverse_strand" if flag & 16 else "flt1_span_forward_strand"] += 1
        hits = set()
        for start, stop, left, right in introns:
            marker = junctions.get((start, stop))
            if marker is None:
                continue
            if min(left, right) >= SENSITIVITY_FLANK:
                flank5[marker] += 1
            if min(left, right) >= PRIMARY_FLANK:
                hits.add(marker)
        if not introns:
            for marker, (start, stop) in bodies.items():
                if pos >= start and end <= stop:
                    hits.add(marker)
        for marker in hits:
            counts[marker] += 1
        if name_hook is not None:
            per_read[fields[0]] = hits
    return counts, flank5, stats, per_read


def hisat2(work, reads, fasta=False, threads=6):
    command = [str(work / "env/bin/hisat2"), "-p", str(threads), "-x", str(work / "ref/grch37/genome"),
               "--known-splicesite-infile", str(work / "ref/splicesites.txt"), "--no-unal"]
    command += ["-f"] if fasta else []
    command += ["-U", str(reads)]
    return subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


parser = argparse.ArgumentParser()
parser.add_argument("--work", required=True, type=Path)
parser.add_argument("--markers", action="store_true")
parser.add_argument("--runs", nargs="*")
args = parser.parse_args()
work = args.work
marker_ids, junctions, bodies = load_markers()
(work / "counts").mkdir(exist_ok=True)

if args.markers:
    process = hisat2(work, work / "ref/marker-50mers.fa", fasta=True)
    counts, flank5, stats, per_read = count_sam(process.stdout, marker_ids, junctions, bodies, name_hook=True)
    log = process.stderr.read()
    process.wait()
    assert process.returncode == 0, log
    totals = {m: 0 for m in marker_ids}
    correct = {m: 0 for m in marker_ids}
    with (work / "ref/marker-50mers.fa").open() as handle:
        names = [line[1:].strip() for line in handle if line.startswith(">")]
    for name in names:
        marker = name.split("|")[0]
        totals[marker] += 1
        if marker in per_read.get(name, ()):
            correct[marker] += 1
    result = {
        "method": "HISAT2 full GRCh37 genome, known splice sites Ensembl GRCh37.87, default scoring; "
                  "a 50-mer is correct if uniquely aligned (NH=1) and counted for its own marker "
                  "under the primary rule",
        "markers": {m: {"synthetic_50mers": totals[m], "unique_and_correct": correct[m],
                        "fraction": round(correct[m] / totals[m], 4) if totals[m] else None}
                    for m in marker_ids},
    }
    (HERE / "hisat2-marker-mappability.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
else:
    runs = args.runs or sorted(p.name.split(".")[0] for p in (work / "cand").glob("*.fastq.gz"))
    for run in runs:
        out = work / "counts" / f"{run}.json"
        if out.exists():
            continue
        process = hisat2(work, work / "cand" / f"{run}.fastq.gz")
        counts, flank5, stats, _ = count_sam(process.stdout, marker_ids, junctions, bodies)
        log = process.stderr.read()
        process.wait()
        assert process.returncode == 0, log
        stage1 = json.loads((work / "logs" / f"{run}.json").read_text())
        result = {"run": run, "total_reads": stage1["total_reads"], "candidate_reads": stage1["candidate_reads"],
                  "hisat2_unique_counts_flank10": counts, "hisat2_unique_junction_counts_flank5": flank5,
                  "exact_50mer_counts": stage1["exact_marker_counts"], "alignment_stats": stats,
                  "hisat2_summary": log.strip().splitlines()}
        out.write_text(json.dumps(result, indent=2) + "\n")
        print(run, {m: c for m, c in counts.items() if c}, flush=True)
