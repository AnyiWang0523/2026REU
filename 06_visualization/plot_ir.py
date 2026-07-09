# ============================================================
# Monthly Supply + PDSI Line Plot - Irrigation
# Each PNG: one HUC4, 2000-2020 all years
# Layout: 3 subplots (PDSI top, GW middle, SW bottom)
# X-axis: 2000/01 - 2020/12 (252 months)
# Y-axis: PDSI / Supply (Mgal/month)
# ============================================================

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import requests
import time
import os
from datetime import datetime

# ============================================================
# Path Configuration
# ============================================================

BASE_DIR   = r"C:\Users\anyiw\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI"
GW_CSV     = os.path.join(BASE_DIR, "output", "IR_GW", "monthly", "IR_GW_monthly_2000_2020.csv")
SW_CSV     = os.path.join(BASE_DIR, "output", "IR_SW", "monthly", "IR_SW_monthly_2000_2020.csv")
OUTPUT_DIR = os.path.join(BASE_DIR, "output", "IR_plots")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================
# STEP 1: Load GW and SW Data
# ============================================================
print("=" * 50)
print("STEP 1: Loading IR GW and SW monthly data...")
print("=" * 50)

gw = pd.read_csv(GW_CSV, dtype={"HUC4": str})
sw = pd.read_csv(SW_CSV, dtype={"HUC4": str})

gw["HUC4"] = gw["HUC4"].str.zfill(4)
sw["HUC4"] = sw["HUC4"].str.zfill(4)

gw["Date"] = pd.to_datetime(gw[["Year","Month"]].assign(Day=1))
sw["Date"] = pd.to_datetime(sw[["Year","Month"]].assign(Day=1))

all_huc4 = sorted(set(gw["HUC4"].unique()) | set(sw["HUC4"].unique()))
print(f"✓ GW records: {len(gw)}")
print(f"✓ SW records: {len(sw)}")
print(f"✓ HUC4 units: {len(all_huc4)}")

# ============================================================
# STEP 2: Load HUC4 -> State Mapping
# ============================================================
print("\n" + "=" * 50)
print("STEP 2: Loading HUC4-State mapping...")
print("=" * 50)

wsdi_csv = os.path.join(BASE_DIR, "output", "IR_GW", "WSDI_IR_GW_HUC4_2000_2020.csv")
wsdi_df  = pd.read_csv(wsdi_csv, dtype={"HUC4": str})
wsdi_df["HUC4"] = wsdi_df["HUC4"].str.zfill(4)
huc4_state = wsdi_df[["HUC4","STATE"]].drop_duplicates().set_index("HUC4")["STATE"].to_dict()
print(f"✓ HUC4-State mapping: {len(huc4_state)} entries")

# ============================================================
# STEP 3: Fetch PDSI for all states
# ============================================================
print("\n" + "=" * 50)
print("STEP 3: Fetching PDSI from NOAA API...")
print("=" * 50)

NOAA_CODES = {
    "AL":"01","AZ":"02","AR":"03","CA":"04","CO":"05","CT":"06",
    "DE":"07","FL":"08","GA":"09","ID":"10","IL":"11","IN":"12",
    "IA":"13","KS":"14","KY":"15","LA":"16","ME":"17","MD":"18",
    "MA":"19","MI":"20","MN":"21","MS":"22","MO":"23","MT":"24",
    "NE":"25","NV":"26","NH":"27","NJ":"28","NM":"29","NY":"30",
    "NC":"31","ND":"32","OH":"33","OK":"34","OR":"35","PA":"36",
    "RI":"37","SC":"38","SD":"39","TN":"40","TX":"41","UT":"42",
    "VT":"43","VA":"44","WA":"45","WV":"46","WI":"47","WY":"48"
}

def fetch_pdsi_annual(state_abbr):
    """Fetch annual PDSI from NOAA, return dict {year: pdsi_value}"""
    code = NOAA_CODES.get(state_abbr.upper())
    if not code:
        return {}
    url = (
        f"https://www.ncei.noaa.gov/access/monitoring/climate-at-a-glance/"
        f"statewide/time-series/{code}/pdsi/all/0/2000-2020/data.json"
    )
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        return {int(k[:4]): float(v["value"]) for k, v in r.json()["data"].items()}
    except Exception as e:
        print(f"  ✗ {state_abbr}: {e}")
        return {}

states_needed = list(set(huc4_state.values()))
pdsi_by_state = {}
for st in states_needed:
    pdsi_by_state[st] = fetch_pdsi_annual(st)
    print(f"  ✓ {st}")
    time.sleep(0.3)

print(f"✓ PDSI fetched for {len(pdsi_by_state)} states")

# ============================================================
# STEP 4: Generate Plots
# ============================================================
print("\n" + "=" * 50)
print("STEP 4: Generating plots...")
print("=" * 50)

all_dates = pd.date_range("2000-01", "2020-12", freq="MS")

total = len(all_huc4)
count = 0

for huc4 in all_huc4:

    huc4_dir = os.path.join(OUTPUT_DIR, huc4)
    os.makedirs(huc4_dir, exist_ok=True)

    state = huc4_state.get(huc4, None)

    # PDSI: expand annual value to monthly
    pdsi_monthly = []
    if state and state in pdsi_by_state:
        for yr in range(2000, 2021):
            val = pdsi_by_state[state].get(yr, None)
            for mo in range(1, 13):
                pdsi_monthly.append({
                    "Date": datetime(yr, mo, 1),
                    "PDSI": val
                })
    pdsi_df = pd.DataFrame(pdsi_monthly)
    if not pdsi_df.empty:
        pdsi_df["Date"] = pd.to_datetime(pdsi_df["Date"])

    # GW and SW for this HUC4
    gw_sub = gw[gw["HUC4"] == huc4].sort_values("Date")
    sw_sub = sw[sw["HUC4"] == huc4].sort_values("Date")

    gw_full = gw_sub.set_index("Date").reindex(all_dates)["GW_Supply_Mgal"].fillna(0)
    sw_full = sw_sub.set_index("Date").reindex(all_dates)["SW_Supply_Mgal"].fillna(0)

    # Plot
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1,
        figsize=(16, 10),
        sharex=True,
        gridspec_kw={"hspace": 0.4}
    )

    fig.suptitle(
        f"Irrigation Monthly Water Use — HUC4: {huc4} | State: {state} | 2000-2020",
        fontsize=13, fontweight="bold", y=0.99
    )

    # PDSI subplot (top)
    if not pdsi_df.empty:
        ax1.plot(pdsi_df["Date"], pdsi_df["PDSI"],
                 color="#4dac26", linewidth=1.5, label="PDSI")
        ax1.axhline(y=0, color="black", linewidth=0.8, linestyle="--", alpha=0.6)
        ax1.fill_between(pdsi_df["Date"], pdsi_df["PDSI"], 0,
                         where=(pdsi_df["PDSI"] >= 0),
                         alpha=0.15, color="#4dac26", label="Wet")
        ax1.fill_between(pdsi_df["Date"], pdsi_df["PDSI"], 0,
                         where=(pdsi_df["PDSI"] < 0),
                         alpha=0.15, color="#d73027", label="Dry")
    ax1.set_title("PDSI (Palmer Drought Severity Index)", fontsize=10)
    ax1.set_ylabel("PDSI", fontsize=9)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(True, linestyle="--", alpha=0.4)

    # GW subplot (middle)
    ax2.plot(all_dates, gw_full.values, color="#2166ac",
             linewidth=1.5, label="Groundwater")
    ax2.fill_between(all_dates, gw_full.values, alpha=0.1, color="#2166ac")
    ax2.set_title("Groundwater Supply (GW)", fontsize=10)
    ax2.set_ylabel("Supply (Mgal/month)", fontsize=9)
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax2.legend(fontsize=8, loc="upper right")
    ax2.grid(True, linestyle="--", alpha=0.4)

    # SW subplot (bottom)
    ax3.plot(all_dates, sw_full.values, color="#d6604d",
             linewidth=1.5, label="Surface Water")
    ax3.fill_between(all_dates, sw_full.values, alpha=0.1, color="#d6604d")
    ax3.set_title("Surface Water Supply (SW)", fontsize=10)
    ax3.set_ylabel("Supply (Mgal/month)", fontsize=9)
    ax3.set_xlabel("Year", fontsize=9)
    ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax3.xaxis.set_major_locator(mdates.YearLocator())
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax3.legend(fontsize=8, loc="upper right")
    ax3.grid(True, linestyle="--", alpha=0.4)

    # Save
    out_png = os.path.join(huc4_dir, f"IR_{huc4}_2000_2020.png")
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)

    count += 1
    if count % 20 == 0:
        print(f"  Progress: {count}/{total} HUC4s")

print(f"\n✓ All {total} plots generated!")
print(f"  Saved in: {OUTPUT_DIR}")
print("\n" + "=" * 50)
print("All done! Results are in output/IR_plots/ folder")
print("=" * 50)