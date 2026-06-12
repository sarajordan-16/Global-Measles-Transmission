#!/usr/bin/env python3
"""
Bar chart of sequences per sampling year, color-coded by continent.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import country_converter as coco

# ── Load metadata ──────────────────────────────────────────────────────────────
df = pd.read_csv(
    "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/metadata_global_measles_data.tsv",
    sep="\t",
    usecols=["accessionVersion", "sampleCollectionDate", "geoLocCountry"],
    dtype=str,
)

# ── Keep only rows with a country ─────────────────────────────────────────────
df = df[df["geoLocCountry"].notna() & (df["geoLocCountry"].str.strip() != "")]

# ── Convert country to continent ───────────────────────────────────────────────
cc = coco.CountryConverter()

def country_to_continent(country_name):
    """Convert a country name to a continent label, separating North and South America."""
    if pd.isna(country_name) or str(country_name).strip() == "":
        return None

    country_name = str(country_name).strip()
    special_cases = {
        "Hong Kong": "Asia",
        "Viet Nam": "Asia",
        "American Samoa": "Oceania",
        "New Caledonia": "Oceania",
        "Gibraltar": "Europe",
        "Serbia and Montenegro": "Europe",
    }
    if country_name in special_cases:
        return special_cases[country_name]

    continent = cc.convert(names=country_name, to="continent")

    if continent == "not found":
        continent = special_cases.get(country_name, "Unknown")
    elif continent in ["America", "Americas"]:
        # Separate North and South America
        north_america_countries = {
            "Canada", "USA", "United States", "Mexico",
            "Belize", "Costa Rica", "El Salvador", "Guatemala",
            "Honduras", "Nicaragua", "Panama",
            "Bahamas", "Barbados", "Bermuda", "Cuba",
            "Dominica", "Dominican Republic", "Grenada",
            "Haiti", "Jamaica", "Saint Kitts and Nevis",
            "Saint Lucia", "Saint Vincent and the Grenadines",
            "Trinidad and Tobago", "Antigua and Barbuda"
        }
        south_america_countries = {
            "Argentina", "Bolivia", "Brazil", "Chile", "Colombia",
            "Ecuador", "Guyana", "Paraguay", "Peru", "Suriname",
            "Uruguay", "Venezuela"
        }
        if country_name in north_america_countries:
            continent = "North America"
        elif country_name in south_america_countries:
            continent = "South America"
        else:
            # Default to North America for unknown Americas countries
            continent = "North America"

    return continent

# ── Extract year ───────────────────────────────────────────────────────────────
def extract_year(date_str):
    """Return 4-digit year string if parseable, else None."""
    if pd.isna(date_str) or str(date_str).strip() == "":
        return None
    parts = str(date_str).strip().split("-")
    year = parts[0]
    if len(year) == 4 and year.isdigit():
        return int(year)
    return None

df["year"] = df["sampleCollectionDate"].apply(extract_year)
df = df[df["year"].notna()].copy()
df["year"] = df["year"].astype(int)
df["continent"] = df["geoLocCountry"].apply(country_to_continent)
df = df[df["continent"].notna() & (df["continent"] != "Unknown")].copy()

# ── Filter for years >= 1990 ──────────────────────────────────────────────────
df = df[df["year"] >= 1990].copy()

# ── Pivot: rows = year, columns = continent ───────────────────────────────────
pivot = (
    df.groupby(["year", "continent"])
    .size()
    .unstack(fill_value=0)
)

# Sort years
pivot = pivot.sort_index()

# ── Assign colours ────────────────────────────────────────────────────────────
continents = list(pivot.columns)
n = len(continents)

# Use a combination of colormaps to get enough distinct colours
cmap1 = plt.colormaps["tab20"].resampled(20)
cmap2 = plt.colormaps["tab20b"].resampled(20)
cmap3 = plt.colormaps["tab20c"].resampled(20)

all_colors = (
    [cmap1(i) for i in range(20)]
    + [cmap2(i) for i in range(20)]
    + [cmap3(i) for i in range(20)]
)
color_map = {continent: all_colors[i % len(all_colors)] for i, continent in enumerate(continents)}

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(20, 9))

years = pivot.index.tolist()
x = np.arange(len(years))
bar_width = 0.8

bottoms = np.zeros(len(years))
for continent in continents:
    vals = pivot[continent].values.astype(float)
    ax.bar(
        x,
        vals,
        bar_width,
        bottom=bottoms,
        color=color_map[continent],
        label=continent,
        linewidth=0,
    )
    bottoms += vals

# ── Axes formatting ───────────────────────────────────────────────────────────
ax.set_xticks(x)
ax.set_xticklabels(years, rotation=45, ha="right", fontsize=9)
ax.set_xlabel("Sampling Year", fontsize=13, labelpad=10)
ax.set_ylabel("Number of Sequences", fontsize=13, labelpad=10)
ax.set_title("All Available Global Measles Sequences by Sampling Year and Continent", fontsize=15, fontweight="bold", pad=15)
ax.yaxis.grid(True, linestyle="--", alpha=0.5)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# ── Legend ────────────────────────────────────────────────────────────────────
# Only include continents that actually appear (total > 0)
active_continents = [c for c in continents if pivot[c].sum() > 0]
handles = [mpatches.Patch(color=color_map[c], label=c) for c in active_continents]

legend = ax.legend(
    handles=handles,
    title="Continent",
    bbox_to_anchor=(1.01, 1),
    loc="upper left",
    fontsize=12,
    title_fontsize=14,
    ncol=2,
    frameon=True,
    borderpad=0.8,
)

plt.tight_layout(rect=[0, 0, 0.72, 1])

out_path = "/Users/sarajordan/Desktop/Gobal_Measles_DTA_analysis/sequences_by_year_continent.png"
plt.savefig(out_path, dpi=200, bbox_inches="tight")
plt.savefig(out_path.replace(".png", ".pdf"), bbox_inches="tight")
print(f"Saved: {out_path}")
print(f"Years covered : {years[0]} – {years[-1]}")
print(f"Total continents : {len(active_continents)}")
print(f"Total sequences plotted: {int(bottoms.sum())}")
