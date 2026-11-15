"""Export the pre-specified FLT1 marker 50-mers and the FLT1 cDNA bait for the compiled stream filter.

The marker 50-mers are produced by the *same* function as the group-blinded pilot
(`scan_streamed_fastq_markers.build_marker_kmers`), so the compiled filter counts exactly the
k-mers the pilot counted. The bait is the cDNA of the five annotated GRCh37 FLT1 transcripts.

Usage:
  python export_marker_kmers.py --cdna-fasta-gz Homo_sapiens.GRCh37.cdna.all.fa.gz --outdir WORK/ref
"""

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Load the pilot scanner's functions without executing its command-line section.
_source = (HERE / "scan_streamed_fastq_markers.py").read_text().split("\nparser = argparse")[0]
_ns = {"__file__": str(HERE / "scan_streamed_fastq_markers.py"), "__name__": "pilot_scanner"}
exec(compile(_source, "scan_streamed_fastq_markers.py", "exec"), _ns)

parser = argparse.ArgumentParser()
parser.add_argument("--cdna-fasta-gz", required=True, type=Path)
parser.add_argument("--outdir", required=True, type=Path)
args = parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)

lookup, unique_by_marker = _ns["build_marker_kmers"](args.cdna_fasta_gz)
with (args.outdir / "marker-kmers.tsv").open("w") as handle:
    for kmer in sorted(lookup):
        for marker in sorted(lookup[kmer]):
            handle.write(f"{marker}\t{kmer}\n")

# Synthetic marker reads (one orientation each) for the aligner mappability confirmation.
with (args.outdir / "marker-50mers.fa").open("w") as handle:
    for marker, kmers in sorted(unique_by_marker.items()):
        for index, kmer in enumerate(sorted(kmers)):
            handle.write(f">{marker}|{index}\n{kmer}\n")

gene = json.loads((HERE / "grch37-flt1-gene.json").read_text())
transcripts = {tx["id"] for tx in gene["Transcript"]}
written = 0
with (args.outdir / "flt1-bait.fa").open("w") as handle:
    for header, sequence in _ns["fasta_records"](args.cdna_fasta_gz):
        transcript = header.split()[0].split(".")[0]
        if transcript in transcripts:
            handle.write(f">{transcript}\n{sequence}\n")
            written += 1
assert written == len(transcripts) == 5

print(json.dumps({m: len(k) for m, k in sorted(unique_by_marker.items())}, indent=2))
