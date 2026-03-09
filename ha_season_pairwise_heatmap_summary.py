from Bio import SeqIO
import pandas as pd
from datetime import datetime
import re
from collections import Counter
from tqdm import tqdm  # for progress bar

# ------------------------
# 1. Load metadata
# ------------------------
meta = pd.read_csv("metadata_JapanKorea.csv", dtype=str)
meta.columns = [c.lower().strip() for c in meta.columns]

# Clean up key column
meta["isolate_id"] = meta["isolate_id"].astype(str).str.strip()
meta = meta.set_index("isolate_id")

# ------------------------
# 2. Helper: normalize country names
# ------------------------
def normalize_country(location):
    """Return 'japan' or 'korea' based on the location string."""
    if pd.isna(location):
        return "unknown"
    loc = location.lower()
    if any(k in loc for k in ["japan", "jp"]):
        return "japan"
    if any(k in loc for k in ["korea", "south korea", "republic of korea", "kr"]):
        return "korea"
    return "unknown"

# ------------------------
# 3. Helper: define flu season
# ------------------------
def get_season(date_str):
    """Return flu season in YY_YY format."""
    if pd.isna(date_str) or not re.match(r"\d{4}", str(date_str)):
        return "Unknown"
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
        try:
            date = datetime.strptime(date_str, fmt)
            break
        except ValueError:
            continue
    else:
        return "Unknown"

    year, month = date.year, getattr(date, "month", 1)
    if month >= 9:
        return f"{str(year)[-2:]}_{str(year+1)[-2:]}"
    else:
        return f"{str(year-1)[-2:]}_{str(year)[-2:]}"

# ------------------------
# 4. Append metadata info to FASTA headers
# ------------------------
input_fasta = "HA_aligned_edit.fasta"
output_fasta = "HA_jpnkor_aligned_edit_appended.fasta"

records_out = []
not_found = []

print("🔍 Matching FASTA records with metadata...")
for rec in tqdm(list(SeqIO.parse(input_fasta, "fasta"))):
    original_id = rec.id.strip()
    match_id = None

    # exact match first
    if original_id in meta.index:
        match_id = original_id
    else:
        # relaxed substring match (avoid partial overlap errors)
        for key in meta.index:
            if re.search(rf"\b{re.escape(key)}\b", original_id):
                match_id = key
                break

    if match_id:
        location = meta.loc[match_id, "location"] if "location" in meta.columns else "Unknown"
        country = normalize_country(location)
        date = meta.loc[match_id, "collection_date"] if "collection_date" in meta.columns else "Unknown"
        season = get_season(date)

        rec.id = f"{original_id}|{country}|{date}|{season}"
        rec.description = ""
        records_out.append(rec)
    else:
        not_found.append(original_id)

SeqIO.write(records_out, output_fasta, "fasta")
print(f"\n✅ Renamed FASTA written: {output_fasta}")
print(f"ℹ️ {len(records_out)} records updated, {len(not_found)} skipped (no metadata match).")

# Log unmatched IDs for inspection
if not_found:
    pd.Series(not_found, name="unmatched_fasta_ids").to_csv("unmatched_ids.csv", index=False)
    print(f"⚠️ Saved {len(not_found)} unmatched IDs to unmatched_ids.csv")

# ------------------------
# 5. Count sequences per country × season
# ------------------------
counts = Counter()
for r in records_out:
    parts = r.id.split("|")
    if len(parts) >= 4:
        country, season = parts[-3], parts[-1]
        counts[(country, season)] += 1

print("\n=== Sequence counts per country × season ===")
for (country, season), n in sorted(counts.items()):
    print(f"{country}\t{season}\t{n}")

# Optional: save counts to file
pd.DataFrame(
    [(*k, v) for k, v in counts.items()],
    columns=["country", "season", "count"]
).to_csv("HA_jpnkor_seq_counts.csv", index=False)
print("\n✅ Saved count summary as HA_jpnkor_seq_counts.csv")
