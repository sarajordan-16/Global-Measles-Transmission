#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
import re
from pathlib import Path

import baltic as bt


DEFAULT_TREE = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/Subsampled_DTA_Analysis/SUBSAMPLED REINTRODUCTION CLADES/combined_with_trait_trees_MCC.tree")
DEFAULT_OUTPUT = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/reintroduction_clades/persistence_times_mcc.csv")
DEFAULT_SUMMARY_OUTPUT = Path("/Users/sarajordan/Desktop/Global_DTA_Targeted/reintroduction_clades/persistence_times_mcc_by_continent.csv")
TRAIT_NAME = "location"
TREE_STRING_REGEX = r"tree TREE_MCC_median"


def date_to_decimal(date_str: str) -> float | None:
    try:
        parsed = dt.datetime.strptime(date_str, "%Y-%m-%d")
    except ValueError:
        return None

    start = dt.datetime(parsed.year, 1, 1)
    end = dt.datetime(parsed.year + 1, 1, 1)
    return parsed.year + (parsed - start).days / (end - start).days


def decimal_to_datestr(decimal_year: float) -> str:
    year = int(decimal_year)
    fraction = decimal_year - year
    days_in_year = 366 if (year % 4 == 0 and (year % 100 != 0 or year % 400 == 0)) else 365
    day_offset = min(round(fraction * days_in_year), days_in_year - 1)
    date_value = dt.datetime(year, 1, 1) + dt.timedelta(days=day_offset)
    return date_value.strftime("%Y-%m-%d")


def extract_tip_date_from_label(label: str) -> str | None:
    """Extract YYYY-MM-DD date from taxon label like PP_0034GH6.1/2008-05-21/Asia/520/D4"""
    match = re.search(r"/(\d{4}-\d{2}-\d{2})/", label)
    return match.group(1) if match else None


def location_from_name(name: str) -> str | None:
    parts = name.split("/")
    if len(parts) >= 3:
        return parts[2]
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Traverse the MCC tree, identify state-change introduction clades, and "
            "export conservative persistence times using subtree.root.absoluteTime."
        )
    )
    parser.add_argument("--tree", type=Path, default=DEFAULT_TREE, help="Path to the MCC tree.")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="CSV output path for persistence-time summaries.",
    )
    parser.add_argument(
        "--summary-output",
        type=Path,
        default=DEFAULT_SUMMARY_OUTPUT,
        help="CSV output path for continent-level persistence summaries.",
    )
    parser.add_argument(
        "--tip-date-cutoff",
        type=float,
        default=1990.0,
        help="Exclude tips sampled before this decimal year. Default: 1990.0",
    )
    parser.add_argument(
        "--min-tips",
        type=int,
        default=3,
        help="Only retain introduction clades with at least this many sampled tips. Default: 3",
    )
    return parser.parse_args()


def load_tree(tree_path: Path) -> bt.tree:
    """Load tree and set node times based on tip dates from taxon labels."""
    tree = bt.loadNexus(
        str(tree_path),
        treestring_regex=TREE_STRING_REGEX,
        absoluteTime=False,  # Don't use internal calibration
    )
    tree.traverse_tree()
    
    # Extract tip dates from taxon labels and set them as absoluteTime
    for node in tree.Objects:
        if node.is_leaf():
            date_str = extract_tip_date_from_label(node.name)
            if date_str:
                node.absoluteTime = date_to_decimal(date_str)
    
    # Propagate times up the tree from dated tips.
    # Branch lengths are stored on the child edge (parent -> child), so:
    # parent_time = child_time - child.length.
    # Use the minimum implied parent time so parent is not younger than any child.
    def set_internal_node_times(node):
        if node.is_leaf():
            return node.absoluteTime
        else:
            implied_parent_times = []
            for child in node.children:
                child_time = set_internal_node_times(child)
                implied_parent_times.append(child_time - child.length)
            node.absoluteTime = min(implied_parent_times)
            return node.absoluteTime
    
    set_internal_node_times(tree.root)
    tree.sortBranches()
    return tree


def prune_old_tips(tree: bt.tree, cutoff: float) -> None:
    tips_to_remove: set[int] = set()

    for obj in tree.Objects:
        if obj.branchType != "leaf":
            continue
        parts = obj.name.split("/")
        if len(parts) < 2:
            continue
        decimal_date = date_to_decimal(parts[1])
        if decimal_date is not None and decimal_date < cutoff:
            tips_to_remove.add(obj.index)

    tree.Objects = [obj for obj in tree.Objects if obj.index not in tips_to_remove]

    for obj in tree.Objects:
        if obj.branchType == "node":
            obj.children = [child for child in obj.children if child.index not in tips_to_remove]

    changed = True
    while changed:
        changed = False
        singleton_indices: set[int] = set()
        for obj in tree.Objects:
            if obj.branchType == "node" and obj != tree.root and len(obj.children) == 1:
                singleton_indices.add(obj.index)
                child = obj.children[0]
                parent = obj.parent
                child.length += obj.length
                child.parent = parent
                parent.children = [child if candidate.index == obj.index else candidate for candidate in parent.children]
                changed = True
        tree.Objects = [obj for obj in tree.Objects if obj.index not in singleton_indices]

    tree.traverse_tree()
    tree.sortBranches()


def patch_tip_locations(tree: bt.tree) -> None:
    for obj in tree.Objects:
        if obj.branchType != "leaf":
            continue
        if TRAIT_NAME not in obj.traits:
            location = location_from_name(obj.name)
            if location:
                obj.traits[TRAIT_NAME] = location


def summarize_introductions(tree: bt.tree, min_tips: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    event_id = 0

    for obj in tree.Objects:
        parent = obj.parent
        if parent is None:
            continue

        destination = obj.traits.get(TRAIT_NAME)
        if destination is None:
            continue

        source = parent.traits.get(TRAIT_NAME, "ancestor")
        if destination == source:
            continue

        traverse_condition = lambda node, state=destination: node.traits.get(TRAIT_NAME) == state
        subtree = tree.subtree(obj, traverse_condition=traverse_condition)
        if subtree is None:
            continue

        subtree.traverse_tree()
        subtree.sortBranches()
        tips = [node for node in subtree.Objects if node.branchType == "leaf"]
        if len(tips) < min_tips:
            continue

        event_id += 1
        intro_time = subtree.root.absoluteTime
        tip_times = [tip.absoluteTime for tip in tips]
        end_time = max(tip_times)
        first_tip_time = min(tip_times)
        persistence_days = max(0.0, (end_time - intro_time) * 365.25)

        rows.append(
            {
                "event_id": event_id,
                "source_location": source,
                "destination_location": destination,
                "introduction_time_decimal": round(intro_time, 6),
                "introduction_time_date": decimal_to_datestr(intro_time),
                "first_descendant_tip_time_decimal": round(first_tip_time, 6),
                "first_descendant_tip_time_date": decimal_to_datestr(first_tip_time),
                "last_descendant_tip_time_decimal": round(end_time, 6),
                "last_descendant_tip_time_date": decimal_to_datestr(end_time),
                "persistence_days": round(persistence_days, 2),
                "tip_count": len(tips),
                "subtree_size": len(subtree.Objects),
            }
        )

    return rows


def highest_posterior_density_interval(values: list[float], mass: float = 0.95) -> tuple[float, float]:
    if not values:
        raise ValueError("Cannot compute HPDI for an empty list.")

    sorted_values = sorted(values)
    sample_size = len(sorted_values)

    if sample_size == 1:
        return sorted_values[0], sorted_values[0]

    window_size = max(1, math.ceil(mass * sample_size))
    if window_size >= sample_size:
        return sorted_values[0], sorted_values[-1]

    best_start = 0
    best_width = float("inf")
    for start in range(0, sample_size - window_size + 1):
        end = start + window_size - 1
        width = sorted_values[end] - sorted_values[start]
        if width < best_width:
            best_width = width
            best_start = start

    best_end = best_start + window_size - 1
    return sorted_values[best_start], sorted_values[best_end]


def summarize_by_continent(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped_months: dict[str, list[float]] = {}
    moves_out_counts: dict[str, int] = {}

    for row in rows:
        continent = str(row["destination_location"])
        source_continent = str(row["source_location"])
        persistence_days = float(row["persistence_days"])
        persistence_months = persistence_days / 30.4375
        grouped_months.setdefault(continent, []).append(persistence_months)
        if source_continent != "ancestor":
            moves_out_counts[source_continent] = moves_out_counts.get(source_continent, 0) + 1

    summary_rows: list[dict[str, object]] = []
    for continent in sorted(grouped_months):
        values = grouped_months[continent]
        median_months = sorted(values)[len(values) // 2] if len(values) % 2 == 1 else (
            sorted(values)[len(values) // 2 - 1] + sorted(values)[len(values) // 2]
        ) / 2
        hpdi_low, hpdi_high = highest_posterior_density_interval(values, mass=0.95)

        summary_rows.append(
            {
                "continent": continent,
                "reintroduction_events": len(values),
                "moves_out_of_continent": moves_out_counts.get(continent, 0),
                "net_movement": len(values) - moves_out_counts.get(continent, 0),
                "median_persistence_months": round(median_months, 2),
                "hpdi_95_low_months": round(hpdi_low, 2),
                "hpdi_95_high_months": round(hpdi_high, 2),
            }
        )

    return summary_rows


def write_csv(rows: list[dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(rows, key=lambda row: float(row["persistence_days"]), reverse=True)
    fieldnames = [
        "event_id",
        "source_location",
        "destination_location",
        "introduction_time_decimal",
        "introduction_time_date",
        "first_descendant_tip_time_decimal",
        "first_descendant_tip_time_date",
        "last_descendant_tip_time_decimal",
        "last_descendant_tip_time_date",
        "persistence_days",
        "tip_count",
        "subtree_size",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted_rows)


def write_summary_csv(rows: list[dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "continent",
        "reintroduction_events",
        "moves_out_of_continent",
        "net_movement",
        "median_persistence_months",
        "hpdi_95_low_months",
        "hpdi_95_high_months",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    tree = load_tree(args.tree)
    prune_old_tips(tree, args.tip_date_cutoff)
    patch_tip_locations(tree)
    rows = summarize_introductions(tree, args.min_tips)
    summary_rows = summarize_by_continent(rows)
    write_csv(rows, args.output)
    write_summary_csv(summary_rows, args.summary_output)
    print(f"Wrote {len(rows)} introduction events to {args.output}")
    print(f"Wrote {len(summary_rows)} continent summaries to {args.summary_output}")


if __name__ == "__main__":
    main()
