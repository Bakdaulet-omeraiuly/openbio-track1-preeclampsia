"""Count pre-specified exact-unique FLT1 50-mers in streamed FASTQ.

The FASTQ must be decompressed and supplied on stdin. This group-blinded pilot is
only a raw-read coverage check; it is not the splice-aware primary analysis.

Example:
  gzip -dc run.fastq.gz | python scan_streamed_fastq_markers.py \
    --cdna-fasta-gz Homo_sapiens.GRCh37.cdna.all.fa.gz --run SRR000000
"""

import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
K = 50
MIN_JUNCTION_FLANK = 10


def fasta_records(path):
    header = None
    chunks = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks).upper()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def canonical_kmer(sequence):
    return min(sequence, reverse_complement(sequence))


def ordered_exons(transcript):
    return sorted(transcript["Exon"], key=lambda exon: exon["start"], reverse=transcript["strand"] < 0)


def build_marker_kmers(cdna_path):
    gene = json.loads((HERE / "grch37-flt1-gene.json").read_text())
    transcript_objects = {tx["id"]: tx for tx in gene["Transcript"]}
    with (HERE / "flt1-marker-definitions.csv").open() as handle:
        markers = list(csv.DictReader(handle))

    flt1_sequences = {}
    for header, sequence in fasta_records(cdna_path):
        transcript = header.split()[0].split(".")[0]
        if transcript in transcript_objects:
            flt1_sequences[transcript] = sequence
    assert set(flt1_sequences) == set(transcript_objects)

    candidates = defaultdict(set)
    for marker in markers:
        transcript = transcript_objects[marker["transcript"]]
        sequence = flt1_sequences[marker["transcript"]]
        exons = ordered_exons(transcript)
        offsets = [0]
        for exon in exons:
            offsets.append(offsets[-1] + exon["end"] - exon["start"] + 1)
        assert offsets[-1] == len(sequence)

        if marker["marker_type"] == "unique_terminal_body_within_flt1":
            terminal = exons[-1]
            start_coordinate = int(marker["start_grch37"])
            end_coordinate = int(marker["end_grch37"])
            if transcript["strand"] < 0:
                local_start = terminal["end"] - end_coordinate
                local_end = terminal["end"] - start_coordinate + 1
            else:
                local_start = start_coordinate - terminal["start"]
                local_end = end_coordinate - terminal["start"] + 1
            body = sequence[offsets[-2] + local_start : offsets[-2] + local_end]
            for start in range(len(body) - K + 1):
                kmer = body[start : start + K]
                if "N" not in kmer:
                    candidates[marker["marker_id"]].add(canonical_kmer(kmer))
        else:
            boundary = None
            for index, (left, right) in enumerate(zip(exons, exons[1:])):
                low, high = sorted([left, right], key=lambda exon: exon["start"])
                if (low["end"] + 1, high["start"] - 1) == (
                    int(marker["start_grch37"]),
                    int(marker["end_grch37"]),
                ):
                    boundary = offsets[index + 1]
                    break
            assert boundary is not None, marker["marker_id"]
            for start in range(boundary - (K - MIN_JUNCTION_FLANK), boundary - MIN_JUNCTION_FLANK + 1):
                kmer = sequence[start : start + K]
                if len(kmer) == K and "N" not in kmer:
                    candidates[marker["marker_id"]].add(canonical_kmer(kmer))

    # Keep only k-mers independently shown to occur exactly once in the full cDNA reference.
    occurrence = defaultdict(int)
    all_candidates = set().union(*candidates.values())
    oriented = {}
    for kmer in all_candidates:
        oriented[kmer] = kmer
        oriented[reverse_complement(kmer)] = kmer
    for _, sequence in fasta_records(cdna_path):
        for start in range(len(sequence) - K + 1):
            key = oriented.get(sequence[start : start + K])
            if key is not None:
                occurrence[key] += 1

    unique_by_marker = {
        marker: {kmer for kmer in kmers if occurrence[kmer] == 1}
        for marker, kmers in candidates.items()
    }
    lookup = defaultdict(set)
    for marker, kmers in unique_by_marker.items():
        for kmer in kmers:
            lookup[kmer].add(marker)
            lookup[reverse_complement(kmer)].add(marker)
    return lookup, unique_by_marker


def scan_fastq(lookup):
    counts = defaultdict(int)
    total_reads = 0
    invalid_records = 0
    while True:
        header = sys.stdin.readline()
        if not header:
            break
        sequence = sys.stdin.readline().rstrip().upper()
        plus = sys.stdin.readline()
        quality = sys.stdin.readline().rstrip()
        if not plus or len(sequence) != len(quality):
            invalid_records += 1
            continue
        total_reads += 1
        hit_markers = set()
        for start in range(len(sequence) - K + 1):
            kmer = sequence[start : start + K]
            if "N" not in kmer:
                hit_markers.update(lookup.get(kmer, ()))
        for marker in hit_markers:
            counts[marker] += 1
    return total_reads, invalid_records, counts


parser = argparse.ArgumentParser()
parser.add_argument("--cdna-fasta-gz", required=True, type=Path)
parser.add_argument("--run", required=True)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

lookup, unique_by_marker = build_marker_kmers(args.cdna_fasta_gz)
total_reads, invalid_records, counts = scan_fastq(lookup)
result = {
    "run": args.run,
    "method": "exact unique 50-mer streamed FASTQ pilot; not splice-aware primary quantification",
    "total_reads": total_reads,
    "invalid_fastq_records": invalid_records,
    "marker_counts": {
        marker: {
            "reads_with_marker": counts[marker],
            "eligible_exact_unique_50mers": len(kmers),
        }
        for marker, kmers in sorted(unique_by_marker.items())
    },
}
payload = json.dumps(result, indent=2) + "\n"
if args.output:
    args.output.write_text(payload)
print(payload, end="")
