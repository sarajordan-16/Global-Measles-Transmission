import datetime
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Patch
import matplotlib.ticker as ticker
import baltic as bt

# ── CONFIG — update these two paths only ──────────────────────────────────────
TREE_PATH  = '/Users/sarajordan/Desktop/Global_DTA_Targeted/Subsampled_DTA_Analysis/SUBSAMPLED REINTRODUCTION CLADES/combined_with_trait_trees_MCC.tree'   # <── change this
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

Continents = ['Africa', 'Asia', 'Europe', 'North_America', 'Oceania', 'South_America']

MIN_TIPS    = 3
Y_SCALE     = 0.015
Y_GAP       = 1
tipSize     = 8
branchWidth = 1.5

continent_colors = {
    'Africa':        'goldenrod',
    'Asia':          'steelblue',
    'Europe':        'indianred',
    'North_America': 'seagreen',
    'Oceania':       'mediumpurple',
    'South_America': 'darkorange',
    'ancestor':      'dimgrey',
}

# ── HELPER ────────────────────────────────────────────────────────────────────
def date_to_decimal(date_str):
    try:
        d     = datetime.datetime.strptime(date_str, '%Y-%m-%d')
        start = datetime.datetime(d.year,     1, 1)
        end   = datetime.datetime(d.year + 1, 1, 1)
        return d.year + (d - start).days / (end - start).days
    except ValueError:
        return None

def location_from_name(name):
    parts = name.split('/')
    return parts[2] if len(parts) >= 3 else None

def c_func(k):
    return continent_colors.get(k.traits.get('location', 'ancestor'), 'grey')

# ── 1. LOAD TREE ──────────────────────────────────────────────────────────────
print("Loading tree...")
ll = bt.loadNexus(
    TREE_PATH,
    tip_regex=r'\/(\d{4}-\d{2}-\d{2})',
    treestring_regex=r'tree TREE_MCC_median',
    absoluteTime=True
)
ll.treeStats()
print(f"Root absolute time: {ll.root.absoluteTime}")

# ── 2. REMOVE PRE-1990 TIPS ───────────────────────────────────────────────────
cutoff = 1990.0
tips_to_remove = set()

for l in ll.Objects:
    if l.branchType == 'leaf':
        parts = l.name.split('/')
        if len(parts) >= 2:
            dec = date_to_decimal(parts[1])
            if dec is not None and dec < cutoff:
                print(f"  Removing pre-1990 tip: {l.name}  →  {dec:.3f}")
                tips_to_remove.add(l.index)

ll.Objects = [l for l in ll.Objects if l.index not in tips_to_remove]
for l in ll.Objects:
    if l.branchType == 'node':
        l.children = [c for c in l.children if c.index not in tips_to_remove]

# Collapse singletons
changed = True
while changed:
    changed = False
    singleton_indices = set()
    for l in ll.Objects:
        if l.branchType == 'node' and l != ll.root and len(l.children) == 1:
            singleton_indices.add(l.index)
            child        = l.children[0]
            parent       = l.parent
            child.length += l.length
            child.parent  = parent
            parent.children = [child if c.index == l.index else c
                                for c in parent.children]
            changed = True
    ll.Objects = [l for l in ll.Objects if l.index not in singleton_indices]

ll.traverse_tree()
ll.sortBranches()
print(f"Tips remaining after pruning: {len([l for l in ll.Objects if l.branchType == 'leaf'])}")

# ── 3. PATCH TIPS MISSING LOCATION ANNOTATION ────────────────────────────────
traitName = 'location'
for l in ll.Objects:
    if l.branchType == 'leaf' and traitName not in l.traits:
        c = location_from_name(l.name)
        if c:
            l.traits[traitName] = c

# ── 4. DISCOVER LOCATION VALUES ──────────────────────────────────────────────
locations = sorted(set(
    l.traits[traitName]
    for l in ll.Objects
    if traitName in l.traits
))
print(f"Locations found: {locations}")

# ── 5. CHECK FOR UNEXPECTED LOCATIONS ────────────────────────────────────────
# Warn if tree contains locations not in expected continent list
unexpected = [loc for loc in locations if loc not in Continents]
if unexpected:
    print(f"WARNING: Unexpected location values found: {unexpected}")
    print("  These will appear as grey in plots. Check your trait annotations.")

# ── 6. STORAGE & SUBTREE EXTRACTION ──────────────────────────────────────────
all_keys      = locations + ['ancestor']
subtype_trees = {c: [] for c in all_keys}

for l in ll.Objects:
    k  = l
    kp = l.parent

    kloc = k.traits.get(traitName)
    if kloc is None:
        continue

    kploc = kp.traits.get(traitName, 'ancestor')

    if kloc != kploc:
        traverse_condition = lambda w, kc=kloc: w.traits.get(traitName) == kc
        subtree = ll.subtree(k, traverse_condition=traverse_condition)

        if subtree:
            n_tips = len([x for x in subtree.Objects if x.branchType == 'leaf'])
            if n_tips >= MIN_TIPS:
                subtree.traverse_tree()
                subtree.sortBranches()
                subtype_trees[kloc].append((kploc, subtree))

for c in locations:
    print(f"  {c}: {len(subtype_trees[c])} introduction events (>= {MIN_TIPS} tips)")

# ── 7. SUMMARY TABLE ─────────────────────────────────────────────────────────
print("\n── Introduction event summary ───────────────────────────────────────")
print(f"{'Destination':<18} {'N events':>8} {'Total tips':>12} {'Median tips':>13}")
for c in Continents:
    events = subtype_trees[c]
    if not events:
        print(f"  {c:<16} {'0':>8} {'-':>12} {'-':>13}")
        continue
    tip_counts = [len([x for x in tr.Objects if x.branchType == 'leaf'])
                  for _, tr in events]
    import statistics
    print(f"  {c:<16} {len(events):>8} {sum(tip_counts):>12} {statistics.median(tip_counts):>13.1f}")
print()

# ── 8. SAVE ONE PDF + PNG PER CONTINENT ──────────────────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

x_attr = lambda k: k.absoluteTime

for continent in Continents:
    if not subtype_trees[continent]:
        print(f"Skipping {continent} — no qualifying subtrees")
        continue

    sorted_trees = sorted(
        subtype_trees[continent],
        key=lambda x: (-x[1].root.absoluteTime, len(x[1].Objects))
    )

    continent_y_units = sum(tr[1].ySpan + Y_GAP for tr in sorted_trees)
    fig_height = max(3, continent_y_units * Y_SCALE)

    fig = plt.figure(figsize=(16, fig_height), facecolor='w')
    gs  = gridspec.GridSpec(1, 1)
    ax  = plt.subplot(gs[0], facecolor='w')

    cumulative_y = 0
    for origin, loc_tree in sorted_trees:
        cy     = cumulative_y
        y_attr = lambda k, _cy=cy: k.y + _cy

        loc_tree.plotTree(ax, x_attr=x_attr, y_attr=y_attr,
                          colour=c_func, width=branchWidth)
        loc_tree.plotPoints(ax, x_attr=x_attr, y_attr=y_attr,
                            size=tipSize, colour=c_func, zorder=100)

        # Large dot at subtree root showing source continent
        oriC = continent_colors.get(origin, 'dimgrey')
        oriX = loc_tree.root.absoluteTime - loc_tree.root.length
        oriY = loc_tree.root.y + cumulative_y
        ax.scatter(oriX, oriY, 120, facecolor=oriC,
                   edgecolor='w', lw=1.2, zorder=200)

        cumulative_y += loc_tree.ySpan + Y_GAP

    # Axes — extend left enough that oldest root dot stays visible in saved figures
    oldest_origin = min(
        loc_tree.root.absoluteTime - loc_tree.root.length
        for _, loc_tree in sorted_trees
    )
    x_min = min(1990, oldest_origin - 3)
    ax.set_xlim(x_min, 2026)
    ax.margins(x=0.01)
    ax.set_ylim(-5, cumulative_y)
    ax.set_yticklabels([])
    ax.tick_params(axis='x', size=4, labelsize=14)
    ax.tick_params(axis='y', size=0)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.grid(axis='x', ls='--', alpha=0.4, lw=0.8)
    [ax.spines[s].set_visible(False) for s in ['top', 'right', 'left']]

    ax.set_title(
        f"{continent}  —  {len(sorted_trees)} introduction events (≥{MIN_TIPS} tips)",
        fontsize=16, fontweight='bold', pad=12,
        color=continent_colors.get(continent, 'black')
    )

    # Legend: only show source continents actually present in this figure
    sources_present = set(origin for origin, _ in sorted_trees)
    legend_elements = [
        Patch(facecolor=continent_colors.get(c, 'grey'), label=c.replace('_', ' '))
        for c in ['ancestor'] + Continents
        if c in sources_present or c == continent
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10,
              framealpha=0.9, title='Source Continent', title_fontsize=11)

    fname   = os.path.join(OUTPUT_DIR, f'measles_reintroductions_{continent}.pdf')
    pngname = os.path.join(OUTPUT_DIR, f'measles_reintroductions_{continent}.png')
    plt.tight_layout()
    plt.savefig(fname,   bbox_inches='tight', pad_inches=0.25)
    plt.savefig(pngname, dpi=300, bbox_inches='tight', pad_inches=0.25)
    plt.close()
    print(f"Saved: {fname}")

print("\nAll continents complete.")