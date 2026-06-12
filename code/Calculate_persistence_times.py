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
OUTPUT_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = OUTPUT_DIR / "persistence_times_mcc.csv"
DEFAULT_SUMMARY_OUTPUT = Path(
    "/Users/sarajordan/Desktop/Global_DTA_Targeted/"
    "reintroduction_clades/persistence_times_mcc_by_continent.csv")

TRAIT_NAME        = "location"
TREE_STRING_REGEX = r"tree TREE_MCC_median"
TARGET_CONTINENT  = "North_America"   # exact string as it appears in your annotations
 
# cluster size bins matching your figure
def size_category(n: int) -> str:
    if n == 1:
        return "Singleton"
    elif n <= 4:
        return "2-4"
    elif n <= 9:
        return "5-9"
    else:
        return "10+"
    
def date_to_decimal(date_str: str) -> float:
    """Convert YYYY-MM-DD to decimal year."""
    parsed = dt.datetime.strptime(date_str, "%Y-%m-%d")
    year_start = dt.datetime(parsed.year, 1, 1)
    year_end = dt.datetime(parsed.year + 1, 1, 1)
    days_in_year = (year_end - year_start).days
    return parsed.year + (parsed - year_start).days / days_in_year

def decimal_year_to_date(decimal: float) -> str:
    """Convert a decimal year (e.g. 2023.456) to a YYYY-MM-DD string."""
    year = int(decimal)
    start = dt.datetime(year, 1, 1)
    end   = dt.datetime(year + 1, 1, 1)
    days  = (decimal - year) * (end - start).days
    date  = start + dt.timedelta(days=days)
    return date.strftime("%Y-%m-%d")

def decimal_year_to_month(decimal: float) -> str:
    """Return abbreviated month name for a decimal year."""
    year  = int(decimal)
    start = dt.datetime(year, 1, 1)
    end   = dt.datetime(year + 1, 1, 1)
    days  = (decimal - year) * (end - start).days
    date  = start + dt.timedelta(days=days)
    return date.strftime("%b")   # e.g. "Apr"

def extract_tip_date_from_label(label: str) -> str | None:
    """Extract YYYY-MM-DD date from taxon label like PP_0034GH6.1/2008-05-21/Asia/520/D4"""
    match = re.search(r"/(\d{4}-\d{2}-\d{2})/", label)
    return match.group(1) if match else None

def get_location(node) -> str | None:
    """Return the location trait annotation for a node, or None."""
    traits = getattr(node, "traits", {})
    return traits.get(TRAIT_NAME, None)
 
 
def is_na(node) -> bool:
    return get_location(node) == TARGET_CONTINENT
 
def load_tree(tree_path: Path) -> bt.tree:
    """Load tree and set node times based on tip dates from taxon labels."""
    myTree = bt.loadNexus(
        str(tree_path),
        treestring_regex=TREE_STRING_REGEX,
        absoluteTime=False,  # Don't use internal calibration
    )
    myTree.traverse_tree()
    
    # Extract tip dates from taxon labels and set them as absoluteTime
    tip_dates: dict[int, float] = {}  # node index -> decimal year
    
    for node in myTree.Objects:
        if node.is_leaf():
            date_str = extract_tip_date_from_label(node.name)
            if date_str:
                tip_dates[node.index] = date_to_decimal(date_str)
                node.absoluteTime = tip_dates[node.index]
    
    # Propagate times up the tree from leaf dates.
    # In baltic/Newick, child.length is the branch length from parent -> child.
    # Therefore each child implies: parent_time = child_time - child.length.
    # Use the minimum implied parent_time so the parent is not younger than any child.
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
    
    set_internal_node_times(myTree.root)
    return myTree


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

 # ── Core traversal ─────────────────────────────────────────────────────────────
def find_reintroduction_roots(tree: bt.tree) -> list:
    """
    Return every node that is:
      - annotated as North America
      - whose parent is annotated as something OTHER than North America
        (i.e. the branch leading into this node crosses a continent boundary)
 
    Each such node is the root of one independent reintroduction clade.
    Nested introductions are counted independently: if a sub-clade later
    transitions out of NA and back in, that inner transition also appears
    as a separate root.
    """
    intro_roots = []
 
    for node in tree.Objects:
        if node.is_leaf():
            continue  # leaves can't be clade roots
 
        if not is_na(node):
            continue
 
        parent = node.parent
        if parent is None:
            continue  # skip tree root
 
        parent_loc = get_location(parent)
        if parent_loc is None:
            continue
 
        if parent_loc != TARGET_CONTINENT:
            intro_roots.append(node)
 
    return intro_roots
 
 
def collect_na_tips(clade_root, all_intro_roots_set: set) -> list:
    """
    Collect all NA tips that are descendants of clade_root, but STOP
    descending into any branch that:
      (a) transitions away from North America, OR
      (b) is itself the root of a nested reintroduction clade
          (those are counted independently).
 
    Returns list of tip nodes.
    """
    na_tips = []
    stack   = list(clade_root.children)
 
    while stack:
        node = stack.pop()
 
        if not is_na(node):
            # This branch has left NA — do not count descendants
            continue
 
        if node in all_intro_roots_set and node is not clade_root:
            # Nested reintroduction — counted independently, stop here
            continue
 
        if node.is_leaf():
            na_tips.append(node)
        else:
            stack.extend(node.children)
 
    return na_tips
 

# ── Main ───────────────────────────────────────────────────────────────────────
def main(tree_path: Path = DEFAULT_TREE,
         output_path: Path = DEFAULT_OUTPUT,
         summary_path: Path = DEFAULT_SUMMARY_OUTPUT) -> None:
 
    print(f"Loading tree from:\n  {tree_path}\n")
    tree = load_tree(tree_path)
    removed = prune_tips_by_year(tree, {1919, 1920})
    if removed:
        print(f"Removed {removed} tips sampled in 1919/1920.\n")
    print(f"Tree loaded. {len([n for n in tree.Objects if n.is_leaf()])} tips total.\n")
 
    intro_roots     = find_reintroduction_roots(tree)
    intro_roots_set = set(intro_roots)
    print(f"Found {len(intro_roots)} North America reintroduction clades.\n")
 
    rows = []
    
    for idx, root in enumerate(intro_roots, start=1):
 
        na_tips = collect_na_tips(root, intro_roots_set)
 
        if not na_tips:
            # Singleton: the clade root itself has no NA descendant tips
            # (can happen if all children immediately leave NA)
            cluster_size    = 1
            persistence_days = 0.0
            most_recent_date = decimal_year_to_date(root.absoluteTime)
        else:
            cluster_size = len(na_tips)   # does NOT include the root node itself
 
            tip_times        = [t.absoluteTime for t in na_tips]
            most_recent_time = max(tip_times)
            persistence_days = max(0.0, (most_recent_time - root.absoluteTime) * 365.25)
            most_recent_date = decimal_year_to_date(most_recent_time)
 
        root_date        = decimal_year_to_date(root.absoluteTime)
        intro_month      = decimal_year_to_month(root.absoluteTime)
        intro_year       = decimal_year_to_date(root.absoluteTime)[:4]  # "2023" or "2024"
 
        # Source continent = parent node annotation
        source_continent = get_location(root.parent) if root.parent else "Unknown"
 
        rows.append({
            "clade_id"          : idx,
            "source_continent"  : source_continent,
            "root_date"         : root_date,
            "intro_month"       : intro_month,
            "cluster_size"      : cluster_size,       # NA tips only
            "size_category"     : size_category(cluster_size),
            "persistence_days"  : round(persistence_days, 2),
            "most_recent_tip"   : most_recent_date,
            "root_absolute_time": round(root.absoluteTime, 6),
        })
 
    # ── Per-clade CSV ──────────────────────────────────────────────────────────
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "clade_id", "source_continent", "root_date", "intro_month",
        "cluster_size", "size_category", "persistence_days",
        "most_recent_tip", "root_absolute_time",
    ]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Per-clade output written to:\n  {output_path}\n")
    
    # ── Introduction vs local transmission summary ─────────────────────────────
    total_introductions = len(rows)
    total_na_tips       = sum(r["cluster_size"] for r in rows)
    total_cases         = total_introductions + total_na_tips
 
    intro_rate          = total_introductions / total_cases if total_cases else 0
    local_rate          = total_na_tips       / total_cases if total_cases else 0
 
    print("=" * 50)
    print("INTRODUCTION vs LOCAL TRANSMISSION SUMMARY")
    print("=" * 50)
    print(f"  Total reintroduction events (clades) : {total_introductions}")
    print(f"  Total local transmission cases        : {total_na_tips}")
    print(f"  Total cases                           : {total_cases}")
    print(f"  Introduction rate                     : {intro_rate:.1%}")
    print(f"  Local transmission rate               : {local_rate:.1%}")
    print()
 
    # ── Per-source-continent summary ───────────────────────────────────────────
    from collections import defaultdict
    by_continent: dict[str, dict] = defaultdict(lambda: {
        "n_introductions": 0,
        "n_local_cases"  : 0,
    })
    for r in rows:
        src = r["source_continent"]
        by_continent[src]["n_introductions"] += 1
        by_continent[src]["n_local_cases"]   += r["cluster_size"]
 
    summary_rows = []
    for src, counts in sorted(by_continent.items()):
        n_intro  = counts["n_introductions"]
        n_local  = counts["n_local_cases"]
        n_total  = n_intro + n_local
        summary_rows.append({
            "source_continent"  : src,
            "n_introductions"   : n_intro,
            "n_local_cases"     : n_local,
            "total_cases"       : n_total,
            "introduction_rate" : round(n_intro / n_total, 4) if n_total else 0,
            "local_rate"        : round(n_local / n_total, 4) if n_total else 0,
        })
 
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_fieldnames = [
        "source_continent", "n_introductions", "n_local_cases",
        "total_cases", "introduction_rate", "local_rate",
    ]
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Per-continent summary written to:\n  {summary_path}\n")
 
 
if __name__ == "__main__":
    main()