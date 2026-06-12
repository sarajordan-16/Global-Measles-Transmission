from __future__ import annotations

import csv
import datetime as dt
import re
from pathlib import Path

import baltic as bt


DEFAULT_TREE = Path(
    "/Users/sarajordan/Desktop/Global_DTA_Targeted/"
    "Subsampled_DTA_Analysis/SUBSAMPLED REINTRODUCTION CLADES/"
    "combined_with_trait_trees_MCC.tree"
)
OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
CLUSTER_OUTPUT = OUTPUT_DIR / "north_america_cluster_membership_by_cluster.csv"
LONG_OUTPUT = OUTPUT_DIR / "north_america_cluster_membership_long.csv"
WIDE_OUTPUT = OUTPUT_DIR / "north_america_cluster_membership_wide.csv"
SUMMARY_OUTPUT = OUTPUT_DIR / "north_america_cluster_summary.csv"

TRAIT_NAME = "location"
TREE_STRING_REGEX = r"tree TREE_MCC_median"
TARGET_CONTINENT = "North_America"
YEARS_TO_REMOVE = {1919, 1920}


def date_to_decimal(date_str: str) -> float:
    parsed = dt.datetime.strptime(date_str, "%Y-%m-%d")
    year_start = dt.datetime(parsed.year, 1, 1)
    year_end = dt.datetime(parsed.year + 1, 1, 1)
    return parsed.year + (parsed - year_start).days / (year_end - year_start).days


def decimal_year_to_date(decimal: float) -> str:
    year = int(decimal)
    start = dt.datetime(year, 1, 1)
    end = dt.datetime(year + 1, 1, 1)
    days = (decimal - year) * (end - start).days
    date = start + dt.timedelta(days=days)
    return date.strftime("%Y-%m-%d")


def extract_tip_date_from_label(label: str) -> str | None:
    match = re.search(r"/(\d{4}-\d{2}-\d{2})/", label)
    return match.group(1) if match else None


def infer_tip_location_from_label(label: str) -> str | None:
    parts = label.split("/")
    return parts[2] if len(parts) >= 3 else None


def get_location(node) -> str | None:
    return getattr(node, "traits", {}).get(TRAIT_NAME)


def is_target_location(node) -> bool:
    return get_location(node) == TARGET_CONTINENT


def load_tree(tree_path: Path) -> bt.tree:
    tree = bt.loadNexus(
        str(tree_path),
        treestring_regex=TREE_STRING_REGEX,
        absoluteTime=False,
    )
    tree.traverse_tree()

    for node in tree.Objects:
        if not node.is_leaf():
            continue
        if TRAIT_NAME not in node.traits:
            inferred_location = infer_tip_location_from_label(node.name)
            if inferred_location:
                node.traits[TRAIT_NAME] = inferred_location

        date_str = extract_tip_date_from_label(node.name)
        if date_str:
            node.absoluteTime = date_to_decimal(date_str)

    def set_internal_node_times(node):
        if node.is_leaf():
            return node.absoluteTime

        implied_parent_times = []
        for child in node.children:
            child_time = set_internal_node_times(child)
            implied_parent_times.append(child_time - child.length)
        node.absoluteTime = min(implied_parent_times)
        return node.absoluteTime

    set_internal_node_times(tree.root)
    tree.sortBranches()
    return tree


def prune_tips_by_year(tree: bt.tree, years_to_remove: set[int]) -> int:
    tip_indices_to_remove: set[int] = set()

    for obj in tree.Objects:
        if obj.branchType != "leaf":
            continue
        match = re.search(r"/(\d{4})-\d{2}-\d{2}", obj.name)
        if not match:
            continue
        sample_year = int(match.group(1))
        if sample_year in years_to_remove:
            tip_indices_to_remove.add(obj.index)

    if not tip_indices_to_remove:
        return 0

    tree.Objects = [obj for obj in tree.Objects if obj.index not in tip_indices_to_remove]

    for obj in tree.Objects:
        if obj.branchType == "node":
            obj.children = [child for child in obj.children if child.index not in tip_indices_to_remove]

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
    return len(tip_indices_to_remove)


def find_reintroduction_roots(tree: bt.tree) -> list:
    intro_roots = []

    for node in tree.Objects:
        if not is_target_location(node):
            continue

        parent = node.parent
        if parent is None:
            continue

        parent_location = get_location(parent)
        if parent_location is None:
            continue

        if parent_location != TARGET_CONTINENT:
            intro_roots.append(node)

    intro_roots.sort(
        key=lambda node: (
            getattr(node, "absoluteTime", float("inf")),
            0 if node.is_leaf() else -len(getattr(node, "children", [])),
            node.name if node.is_leaf() else f"node_{node.index}",
        )
    )
    return intro_roots


def collect_cluster_tips(clade_root, all_intro_roots_set: set) -> list:
    if clade_root.is_leaf():
        return [clade_root]

    collected_tips = []
    stack = list(clade_root.children)

    while stack:
        node = stack.pop()

        if not is_target_location(node):
            continue

        if node in all_intro_roots_set and node is not clade_root:
            continue

        if node.is_leaf():
            collected_tips.append(node)
        else:
            stack.extend(node.children)

    collected_tips.sort(key=lambda tip: (tip.absoluteTime, tip.name))
    return collected_tips


def write_long_csv(rows: list[dict], output_path: Path) -> None:
    fieldnames = [
        "cluster_id",
        "source_continent",
        "cluster_root_date",
        "cluster_root_absolute_time",
        "cluster_size",
        "tip_name",
        "tip_date",
        "tip_absolute_time",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_cluster_csv(rows: list[dict], output_path: Path) -> None:
    fieldnames = [
        "cluster_id",
        "source_continent",
        "cluster_root_date",
        "cluster_root_absolute_time",
        "cluster_size",
        "first_tip_date",
        "last_tip_date",
        "tip_names",
        "tip_dates",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_wide_csv(cluster_to_tips: dict[str, list[str]], output_path: Path) -> None:
    cluster_ids = list(cluster_to_tips)
    max_rows = max((len(tips) for tips in cluster_to_tips.values()), default=0)

    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(cluster_ids)
        for row_idx in range(max_rows):
            writer.writerow([
                cluster_to_tips[cluster_id][row_idx] if row_idx < len(cluster_to_tips[cluster_id]) else ""
                for cluster_id in cluster_ids
            ])


def write_summary_csv(rows: list[dict], output_path: Path) -> None:
    fieldnames = [
        "cluster_id",
        "source_continent",
        "cluster_root_date",
        "cluster_root_absolute_time",
        "cluster_size",
        "first_tip_date",
        "last_tip_date",
    ]
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(tree_path: Path = DEFAULT_TREE) -> None:
    print(f"Loading tree from:\n  {tree_path}\n")
    tree = load_tree(tree_path)
    removed = prune_tips_by_year(tree, YEARS_TO_REMOVE)
    if removed:
        print(f"Removed {removed} tips sampled in {sorted(YEARS_TO_REMOVE)}.\n")

    intro_roots = find_reintroduction_roots(tree)
    intro_roots_set = set(intro_roots)
    print(f"Found {len(intro_roots)} North America introduction roots.\n")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    long_rows: list[dict] = []
    cluster_rows: list[dict] = []
    summary_rows: list[dict] = []
    cluster_to_tips: dict[str, list[str]] = {}

    retained_cluster_count = 0

    for root in intro_roots:
        tips = collect_cluster_tips(root, intro_roots_set)

        if not tips:
            continue

        retained_cluster_count += 1
        cluster_id = f"cluster_{retained_cluster_count:03d}"

        root_date = decimal_year_to_date(root.absoluteTime)
        source_continent = get_location(root.parent) if root.parent else "Unknown"
        cluster_to_tips[cluster_id] = [tip.name for tip in tips]
        first_tip_date = decimal_year_to_date(min(tip.absoluteTime for tip in tips))
        last_tip_date = decimal_year_to_date(max(tip.absoluteTime for tip in tips))

        cluster_rows.append(
            {
                "cluster_id": cluster_id,
                "source_continent": source_continent,
                "cluster_root_date": root_date,
                "cluster_root_absolute_time": round(root.absoluteTime, 6),
                "cluster_size": len(tips),
                "first_tip_date": first_tip_date,
                "last_tip_date": last_tip_date,
                "tip_names": ";".join(tip.name for tip in tips),
                "tip_dates": ";".join((extract_tip_date_from_label(tip.name) or "") for tip in tips),
            }
        )

        summary_rows.append(
            {
                "cluster_id": cluster_id,
                "source_continent": source_continent,
                "cluster_root_date": root_date,
                "cluster_root_absolute_time": round(root.absoluteTime, 6),
                "cluster_size": len(tips),
                "first_tip_date": first_tip_date,
                "last_tip_date": last_tip_date,
            }
        )

        for tip in tips:
            tip_date = extract_tip_date_from_label(tip.name)
            long_rows.append(
                {
                    "cluster_id": cluster_id,
                    "source_continent": source_continent,
                    "cluster_root_date": root_date,
                    "cluster_root_absolute_time": round(root.absoluteTime, 6),
                    "cluster_size": len(tips),
                    "tip_name": tip.name,
                    "tip_date": tip_date or "",
                    "tip_absolute_time": round(tip.absoluteTime, 6),
                }
            )

    write_cluster_csv(cluster_rows, CLUSTER_OUTPUT)
    write_long_csv(long_rows, LONG_OUTPUT)
    write_wide_csv(cluster_to_tips, WIDE_OUTPUT)
    write_summary_csv(summary_rows, SUMMARY_OUTPUT)

    print(f"Saved cluster-level membership table:\n  {CLUSTER_OUTPUT}")
    print(f"Saved long membership table:\n  {LONG_OUTPUT}")
    print(f"Saved wide cluster table:\n  {WIDE_OUTPUT}")
    print(f"Saved cluster summary table:\n  {SUMMARY_OUTPUT}")
    print(f"\nRetained {len(summary_rows)} North America clusters with at least one North American tip.")


if __name__ == "__main__":
    main()
