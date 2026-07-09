import geopandas as gpd
import os

# ============================================================
# Path Configuration
# ============================================================

# Path to the downloaded WBD geodatabase
GDB_PATH = r"C:\Users\wyuoc\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI_GW\data\GBD\WBD_National_GDB.gdb"

# Output folder for the extracted HUC4 shapefile (will be created automatically)
OUT_DIR  = r"C:\Users\wyuoc\OneDrive - University of Illinois - Urbana\GIS\2026REU\WSDI_GW\data\HUC4"

# ============================================================
# Extract HUC4 Layer from GDB
# ============================================================

# Create output folder if it doesn't exist
os.makedirs(OUT_DIR, exist_ok=True)

# Read the WBDHU4 layer from the geodatabase
# This may take 1-2 minutes due to file size
print("Reading HUC4 layer from GDB (this may take 1-2 minutes)...")
huc4 = gpd.read_file(GDB_PATH, layer="WBDHU4")

# Print basic info to verify the data loaded correctly
print(f"Success! Total HUC4 units: {len(huc4)}")
print(f"Available fields: {huc4.columns.tolist()}")

# ============================================================
# Save as Shapefile
# ============================================================

# Define output shapefile path
out_path = os.path.join(OUT_DIR, "WBDHU4.shp")

# Save to shapefile
huc4.to_file(out_path)
print(f"✓ Saved to: {out_path}")