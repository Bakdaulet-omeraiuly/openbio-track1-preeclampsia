"""Build FLT1 isoform feasibility tables without testing clinical-group differences."""

from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TRANSCRIPTS = [
    "ENST00000282397",
    "ENST00000541932",
    "ENST00000539099",
    "ENST00000540678",
    "ENST00000543394",
]
UNIPROT = {
    "ENST00000282397": "P17948-1 / Flt1",
    "ENST00000541932": "P17948-3 / sFlt1-14",
    "ENST00000539099": "P17948-4",
    "ENST00000540678": "P17948-5 / i15",
    "ENST00000543394": "P17948-8",
}


def ordered_exons(transcript):
    return sorted(transcript["Exon"], key=lambda exon: exon["start"], reverse=transcript["strand"] < 0)


def exon_key(exon):
    return exon["seq_region_name"], exon["start"], exon["end"], exon["strand"]


gene = json.loads((HERE / "grch37-flt1-gene.json").read_text())
transcript_objects = {tx["id"]: tx for tx in gene["Transcript"]}
assert set(transcript_objects) == set(TRANSCRIPTS), "GRCh37 FLT1 transcript set changed"

exon_users = Counter()
junction_users = Counter()
structures = {}
for transcript_id, transcript in transcript_objects.items():
    exons = ordered_exons(transcript)
    exon_keys = [exon_key(exon) for exon in exons]
    junctions = []
    for left, right in zip(exons, exons[1:]):
        low, high = sorted([left, right], key=lambda exon: exon["start"])
        junctions.append((low["seq_region_name"], low["end"] + 1, high["start"] - 1, transcript["strand"]))
    for key in set(exon_keys):
        exon_users[key] += 1
    for key in set(junctions):
        junction_users[key] += 1
    structures[transcript_id] = (exons, exon_keys, junctions)

count_parts = []
for group, suffix in [("control", "ControlONLY"), ("IUGR", "IUGRONLY"), ("PE", "PEONLY"), ("PE+IUGR", "PEIUGR")]:
    counts = pd.read_csv(ROOT / "data_external" / f"GSE114691_MasterCount_{suffix}.txt.gz", sep="\t", index_col=0)
    selected = counts.loc[TRANSCRIPTS].T
    selected.index = [f"{group}|{sample}" for sample in selected.index]
    selected["library_total"] = counts.sum(axis=0).values
    count_parts.append(selected)
processed = pd.concat(count_parts)
transcript_counts = processed[TRANSCRIPTS]
cpm = transcript_counts.div(processed.library_total, axis=0) * 1e6
naive_share = transcript_counts.div(transcript_counts.sum(axis=1), axis=0)
five_transcript_sum = transcript_counts.sum(axis=1)

rows = []
marker_rows = []
for transcript_id in TRANSCRIPTS:
    exons, exon_keys, junctions = structures[transcript_id]
    terminal = exons[-1]
    terminal_key = exon_keys[-1]
    other_positions = set()
    for other_id, (other_exons, _, _) in structures.items():
        if other_id == transcript_id:
            continue
        for exon in other_exons:
            other_positions.update(range(exon["start"], exon["end"] + 1))
    terminal_positions = set(range(terminal["start"], terminal["end"] + 1))
    unique_positions = sorted(terminal_positions - other_positions)
    longest_contiguous = 0
    current = 0
    previous = None
    for position in unique_positions:
        current = current + 1 if previous is not None and position == previous + 1 else 1
        longest_contiguous = max(longest_contiguous, current)
        previous = position
    terminal_junction = junctions[-1]
    left_exon, right_exon = exons[-2], exons[-1]
    unique_junctions = [junction for junction in junctions if junction_users[junction] == 1]
    body_intervals = []
    if unique_positions:
        interval_start = previous_position = unique_positions[0]
        for position in unique_positions[1:]:
            if position != previous_position + 1:
                body_intervals.append((interval_start, previous_position))
                interval_start = position
            previous_position = position
        body_intervals.append((interval_start, previous_position))
    for index, (start, end) in enumerate(body_intervals, 1):
        marker_rows.append({
            "marker_id": f"{transcript_id}_terminal_body_{index}",
            "transcript": transcript_id,
            "marker_type": "unique_terminal_body_within_flt1",
            "chromosome": terminal["seq_region_name"],
            "start_grch37": start,
            "end_grch37": end,
            "strand": transcript_objects[transcript_id]["strand"],
            "primary_role": "secondary_cross_check" if transcript_id in {"ENST00000282397", "ENST00000541932"} else "exploratory",
            "note": "Full-genome 50-mer mappability still requires validation",
        })
    for index, junction in enumerate(unique_junctions, 1):
        if transcript_id == "ENST00000541932":
            role = "primary_numerator"
        elif transcript_id == "ENST00000282397":
            role = "primary_denominator_component"
        else:
            role = "exploratory"
        marker_rows.append({
            "marker_id": f"{transcript_id}_unique_junction_{index}",
            "transcript": transcript_id,
            "marker_type": "unique_splice_junction_within_flt1",
            "chromosome": junction[0],
            "start_grch37": junction[1],
            "end_grch37": junction[2],
            "strand": junction[3],
            "primary_role": role,
            "note": "Coordinates are STAR-style intron boundaries; full-genome mappability still requires validation",
        })
    rows.append({
        "transcript": transcript_id,
        "uniprot_exact_match": UNIPROT[transcript_id],
        "exon_count": len(exons),
        "terminal_exon_id": terminal["id"],
        "terminal_exon_coordinates_grch37": f'{terminal["seq_region_name"]}:{terminal["start"]}-{terminal["end"]}',
        "terminal_exon_length_nt": terminal["end"] - terminal["start"] + 1,
        "terminal_unique_body_bases_vs_other_flt1": len(unique_positions),
        "longest_contiguous_unique_terminal_body_nt": longest_contiguous,
        "supports_50nt_body_read_within_flt1": longest_contiguous >= 50,
        "unique_exon_count_by_exact_coordinates": sum(exon_users[key] == 1 for key in exon_keys),
        "unique_junction_count_within_flt1": sum(junction_users[key] == 1 for key in junctions),
        "terminal_junction_left_exon": f'{left_exon["seq_region_name"]}:{left_exon["start"]}-{left_exon["end"]}',
        "terminal_junction_right_exon": f'{right_exon["seq_region_name"]}:{right_exon["start"]}-{right_exon["end"]}',
        "terminal_junction_intron_grch37": f"{terminal_junction[0]}:{terminal_junction[1]}-{terminal_junction[2]}:{terminal_junction[3]}",
        "terminal_junction_unique_within_flt1": junction_users[terminal_junction] == 1,
        "unique_junction_introns_grch37": ";".join(f"{junction[0]}:{junction[1]}-{junction[2]}:{junction[3]}" for junction in unique_junctions),
        "processed_median_raw_count": transcript_counts[transcript_id].median(),
        "processed_spearman_with_sum_five_rows": transcript_counts[transcript_id].corr(five_transcript_sum, method="spearman"),
        "naive_raw_count_share_mean_not_isoform_fraction": naive_share[transcript_id].mean(),
        "naive_raw_count_share_sd_not_isoform_fraction": naive_share[transcript_id].std(),
    })

audit = pd.DataFrame(rows)
audit.to_csv(HERE / "grch37-flt1-exon-junction-audit.csv", index=False)
pd.DataFrame(marker_rows).to_csv(HERE / "flt1-marker-definitions.csv", index=False)
np.log2(cpm + 1).corr(method="spearman").to_csv(HERE / "processed-counts-spearman.csv")

runinfo = pd.read_csv(HERE / "sra-runinfo.csv")
summary = pd.DataFrame([{
    "study": "SRP148556",
    "runs": len(runinfo),
    "library_layout": ";".join(sorted(runinfo.LibraryLayout.unique())),
    "median_read_length_nt": runinfo.avgLength.median(),
    "min_read_length_nt": runinfo.avgLength.min(),
    "max_read_length_nt": runinfo.avgLength.max(),
    "median_spots_per_run": runinfo.spots.median(),
    "total_bases": int(runinfo.bases.sum()),
    "total_sra_size_gib": runinfo.size_MB.sum() / 1024,
}])
summary.to_csv(HERE / "sra-study-summary.csv", index=False)

assert np.allclose(transcript_counts.to_numpy(), np.round(transcript_counts.to_numpy()))
assert audit.loc[audit.transcript == "ENST00000541932", "terminal_junction_unique_within_flt1"].item()
assert not audit.loc[audit.transcript == "ENST00000282397", "terminal_junction_unique_within_flt1"].item()
print(audit.to_string(index=False))
print("\nSRA summary\n", summary.to_string(index=False))
