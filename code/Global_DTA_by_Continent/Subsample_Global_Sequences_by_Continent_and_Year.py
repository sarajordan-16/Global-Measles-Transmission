#!/usr/bin/env python3
"""
subsample_fasta_by_continent_year.py
-----------------------------------

Subsample a measles FASTA file using stratified sampling by continent and year.

This script expects FASTA headers to contain fields separated by '/' where the
second field is a date (YYYY-MM-DD) and the third field is the continent
(for example: '>PP_0034BRV.1/2014-04-11/Asia/456/B3'). It groups sequences into
continent x year strata and performs a two-stage allocation:

    1. Allocate the global target sample size across continents with an optional
         cap per continent (``--max-per-continent``).
    2. For each continent, allocate its continental quota across years (strata)
         with an option to enforce a minimum per non-empty stratum (``--min-per-stratum``).

The sampled sequences are written to a FASTA file and a TSV metadata file that
lists the header, date, year, continent and stratum for each sampled record.

Key options:
    --input            Input FASTA file (default set in the script).
    --output           Output FASTA for the subsample.
    --metadata-output  TSV describing sampled records and their strata.
    --target-n         Maximum number of sequences to keep (default 5000).
    --max-per-continent  Per-continent cap (default 1200).
    --min-per-stratum  Minimum to retain per non-empty continent-year stratum when possible.
    --seed             Random seed for reproducible sampling (default 42).

Notes:
    - The script reads all sequences into memory; very large FASTA files may require modifications.
    - Headers must match the expected format, otherwise a ValueError is raised.

Example:
    python3 subsample_fasta_by_continent_year.py --input input.fasta --output out.fasta --target-n 1000

"""

from __future__ import annotations

import argparse
import csv
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Record:
    header: str
    sequence: str
    date: str
    year: int
    continent: str


DEFAULT_INPUT = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/FASTA files/Continents_aligned_N450_final_continents.fasta")
DEFAULT_OUTPUT = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/FASTA files/Continents_aligned_N5000_subsampled.fasta")
DEFAULT_METADATA = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/FASTA files/Continents_aligned_N5000_subsampled_metadata.tsv")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Subsample a measles FASTA by continent and year using stratified sampling. "
            "Headers are expected to look like '>PP_0034BRV.1/2014-04-11/Asia/456/B3'."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input FASTA file.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output FASTA file.")
    parser.add_argument(
        "--metadata-output",
        type=Path,
        default=DEFAULT_METADATA,
        help="TSV describing sampled records and their strata.",
    )
    parser.add_argument("--target-n", type=int, default=5000, help="Maximum number of sequences to keep.")
    parser.add_argument(
        "--max-per-continent",
        type=int,
        default=1200,
        help=(
            "Maximum number of sequences allowed from any one continent. "
            "Default: 1200"
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducible subsampling.",
    )
    parser.add_argument(
        "--min-per-stratum",
        type=int,
        default=1,
        help="Minimum number to retain per non-empty continent-year stratum when possible.",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_lines: list[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)

    if header is not None:
        records.append((header, "".join(seq_lines)))

    return records


def parse_record(header: str, sequence: str) -> Record:
    parts = header.split("/")
    if len(parts) < 3:
        raise ValueError(f"Header does not have expected fields: {header}")

    date = parts[1]
    continent = parts[2]
    try:
        year = int(date.split("-")[0])
    except ValueError as exc:
        raise ValueError(f"Could not parse year from header: {header}") from exc

    return Record(
        header=header,
        sequence=sequence,
        date=date,
        year=year,
        continent=continent,
    )


def largest_remainder_allocation(
    capacities: dict[object, int],
    target_n: int,
    min_per_stratum: int,
) -> dict[object, int]:
    strata = list(capacities)
    if target_n >= sum(capacities.values()):
        return capacities.copy()

    allocation = {key: 0 for key in strata}

    eligible = [key for key, cap in capacities.items() if cap > 0]
    baseline_total = min_per_stratum * len(eligible)
    if baseline_total <= target_n:
        for key in eligible:
            allocation[key] = min(min_per_stratum, capacities[key])
    else:
        # If the target is too small to guarantee the minimum across all strata,
        # give one sequence to as many strata as possible, starting with largest strata.
        for key, _cap in sorted(capacities.items(), key=lambda item: item[1], reverse=True)[:target_n]:
            allocation[key] = 1
        return allocation

    remaining = target_n - sum(allocation.values())
    if remaining <= 0:
        return allocation

    residual_capacities = {key: capacities[key] - allocation[key] for key in strata}
    residual_total = sum(residual_capacities.values())
    if residual_total == 0:
        return allocation

    quotas: dict[tuple[str, int], float] = {}
    floors: dict[tuple[str, int], int] = {}
    for key, residual in residual_capacities.items():
        quota = (residual / residual_total) * remaining
        quotas[key] = quota
        floors[key] = min(int(quota), residual)

    for key, floor_value in floors.items():
        allocation[key] += floor_value

    leftover = target_n - sum(allocation.values())
    if leftover <= 0:
        return allocation

    remainders = sorted(
        strata,
        key=lambda key: (quotas[key] - floors[key], residual_capacities[key]),
        reverse=True,
    )
    for key in remainders:
        if leftover == 0:
            break
        if allocation[key] < capacities[key]:
            allocation[key] += 1
            leftover -= 1

    return allocation


def capped_continent_allocation(
    continent_capacities: dict[str, int],
    target_n: int,
    max_per_continent: int,
) -> dict[str, int]:
    if target_n >= sum(continent_capacities.values()):
        return continent_capacities.copy()

    effective_capacities = {
        continent: min(count, max_per_continent)
        for continent, count in continent_capacities.items()
    }
    if sum(effective_capacities.values()) < target_n:
        raise ValueError(
            "Target sample size cannot be reached with the chosen --max-per-continent cap. "
            f"Target={target_n}, capped capacity={sum(effective_capacities.values())}."
        )

    allocation = {continent: 0 for continent in continent_capacities}
    remaining_target = target_n
    remaining = effective_capacities.copy()

    while remaining_target > 0:
        total_remaining_capacity = sum(remaining.values())
        if total_remaining_capacity == 0:
            break

        quotas = {
            continent: (capacity / total_remaining_capacity) * remaining_target
            for continent, capacity in remaining.items()
        }
        floors = {
            continent: min(int(quota), remaining[continent])
            for continent, quota in quotas.items()
        }

        assigned_this_round = sum(floors.values())
        for continent, n_keep in floors.items():
            allocation[continent] += n_keep
            remaining[continent] -= n_keep

        remaining_target -= assigned_this_round
        if remaining_target <= 0:
            break

        ranked = sorted(
            remaining,
            key=lambda continent: (quotas[continent] - floors[continent], remaining[continent]),
            reverse=True,
        )
        progress = False
        for continent in ranked:
            if remaining_target == 0:
                break
            if remaining[continent] > 0:
                allocation[continent] += 1
                remaining[continent] -= 1
                remaining_target -= 1
                progress = True

        if not progress:
            break

    return allocation


def write_fasta(records: list[Record], output_path: Path, wrap: int = 80) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for record in records:
            handle.write(f">{record.header}\n")
            sequence = record.sequence
            for start in range(0, len(sequence), wrap):
                handle.write(sequence[start : start + wrap] + "\n")


def write_metadata(records: list[Record], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["header", "date", "year", "continent", "stratum"])
        for record in records:
            writer.writerow([record.header, record.date, record.year, record.continent, f"{record.continent}|{record.year}"])


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    parsed_records = [parse_record(header, sequence) for header, sequence in read_fasta(args.input)]
    buckets: dict[tuple[str, int], list[Record]] = defaultdict(list)
    for record in parsed_records:
        buckets[(record.continent, record.year)].append(record)

    continent_buckets: dict[str, list[Record]] = defaultdict(list)
    for record in parsed_records:
        continent_buckets[record.continent].append(record)

    continent_capacities = {
        continent: len(records)
        for continent, records in continent_buckets.items()
    }
    continent_allocation = capped_continent_allocation(
        continent_capacities,
        args.target_n,
        args.max_per_continent,
    )

    allocation: dict[tuple[str, int], int] = {}
    for continent, continent_target in continent_allocation.items():
        year_capacities = {
            key: len(records)
            for key, records in buckets.items()
            if key[0] == continent
        }
        year_allocation = largest_remainder_allocation(
            year_capacities,
            continent_target,
            args.min_per_stratum,
        )
        allocation.update(year_allocation)

    sampled: list[Record] = []
    for key in sorted(buckets):
        bucket = buckets[key]
        keep_n = allocation[key]
        if keep_n == 0:
            continue
        sampled.extend(rng.sample(bucket, keep_n) if keep_n < len(bucket) else bucket)

    sampled.sort(key=lambda record: (record.continent, record.year, record.date, record.header))

    write_fasta(sampled, args.output)
    write_metadata(sampled, args.metadata_output)

    print(f"Input sequences: {len(parsed_records)}")
    print(f"Target sequences: {args.target_n}")
    print(f"Sampled sequences: {len(sampled)}")
    print(f"Max per continent: {args.max_per_continent}")
    print(f"Strata retained: {sum(1 for n in allocation.values() if n > 0)}")
    print(f"Wrote FASTA: {args.output}")
    print(f"Wrote metadata: {args.metadata_output}")


if __name__ == "__main__":
    main()
