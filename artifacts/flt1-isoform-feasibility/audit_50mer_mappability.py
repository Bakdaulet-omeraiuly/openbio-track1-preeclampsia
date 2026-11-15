"""Audit exact 50-mer uniqueness of FLT1 markers against Ensembl GRCh37 cDNA.

Usage:
  python audit_50mer_mappability.py --cdna-fasta-gz /path/Homo_sapiens.GRCh37.cdna.all.fa.gz

This is an exact-match transcriptome screen, not a substitute for aligner-based genome mappability.
"""

import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path

import pandas as pd

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
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def ordered_exons(transcript):
    return sorted(transcript["Exon"], key=lambda exon: exon["start"], reverse=transcript["strand"] < 0)


parser = argparse.ArgumentParser()
parser.add_argument("--cdna-fasta-gz", required=True, type=Path)
args = parser.parse_args()

gene = json.loads((HERE / "grch37-flt1-gene.json").read_text())
transcript_objects = {tx["id"]: tx for tx in gene["Transcript"]}
markers = pd.read_csv(HERE / "flt1-marker-definitions.csv")

flt1_sequences = {}
for header, sequence in fasta_records(args.cdna_fasta_gz):
    transcript = header.split()[0].split(".")[0]
    if transcript in transcript_objects:
        flt1_sequences[transcript] = sequence
assert set(flt1_sequences) == set(transcript_objects)

marker_kmers = defaultdict(set)
for marker in markers.itertuples(index=False):
    transcript = transcript_objects[marker.transcript]
    sequence = flt1_sequences[marker.transcript]
    exons = ordered_exons(transcript)
    lengths = [exon["end"] - exon["start"] + 1 for exon in exons]
    offsets = [0]
    for length in lengths:
        offsets.append(offsets[-1] + length)
    assert offsets[-1] == len(sequence)

    if marker.marker_type == "unique_terminal_body_within_flt1":
        terminal = exons[-1]
        if transcript["strand"] < 0:
            local_start = terminal["end"] - marker.end_grch37
            local_end = terminal["end"] - marker.start_grch37 + 1
        else:
            local_start = marker.start_grch37 - terminal["start"]
            local_end = marker.end_grch37 - terminal["start"] + 1
        body = sequence[offsets[-2] + local_start : offsets[-2] + local_end]
        for start in range(0, len(body) - K + 1):
            kmer = body[start : start + K]
            if "N" not in kmer:
                marker_kmers[marker.marker_id].add(canonical_kmer(kmer))
    else:
        boundary = None
        for index, (left, right) in enumerate(zip(exons, exons[1:])):
            low, high = sorted([left, right], key=lambda exon: exon["start"])
            intron_start, intron_end = low["end"] + 1, high["start"] - 1
            if intron_start == marker.start_grch37 and intron_end == marker.end_grch37:
                boundary = offsets[index + 1]
                break
        assert boundary is not None, marker.marker_id
        for start in range(boundary - (K - MIN_JUNCTION_FLANK), boundary - MIN_JUNCTION_FLANK + 1):
            kmer = sequence[start : start + K]
            if len(kmer) == K and "N" not in kmer:
                marker_kmers[marker.marker_id].add(canonical_kmer(kmer))

query_kmers = set().union(*marker_kmers.values())
oriented_to_canonical = {}
for kmer in query_kmers:
    oriented_to_canonical[kmer] = kmer
    oriented_to_canonical[reverse_complement(kmer)] = kmer
occurrences = defaultdict(int)
transcript_hits = defaultdict(set)
for header, sequence in fasta_records(args.cdna_fasta_gz):
    transcript = header.split()[0]
    for start in range(0, len(sequence) - K + 1):
        key = oriented_to_canonical.get(sequence[start : start + K])
        if key is not None:
            occurrences[key] += 1
            transcript_hits[key].add(transcript)

rows = []
for marker in markers.itertuples(index=False):
    kmers = sorted(marker_kmers[marker.marker_id])
    unique = [kmer for kmer in kmers if occurrences[kmer] == 1]
    rows.append({
        "marker_id": marker.marker_id,
        "transcript": marker.transcript,
        "marker_type": marker.marker_type,
        "primary_role": marker.primary_role,
        "candidate_exact_50mers": len(kmers),
        "unique_exact_50mers_full_grch37_cdna": len(unique),
        "fraction_unique_exact_50mers": len(unique) / len(kmers) if kmers else 0,
        "max_transcriptome_occurrences": max((occurrences[kmer] for kmer in kmers), default=0),
        "passes_exact_transcriptome_mappability": len(unique) > 0,
        "example_unique_50mer": unique[0] if unique else "",
        "limitation": "Exact 50-mer cDNA screen only; confirm with the selected splice-aware genome aligner",
    })

result = pd.DataFrame(rows)
result.to_csv(HERE / "marker-50mer-mappability.csv", index=False)
print(result.to_string(index=False))
