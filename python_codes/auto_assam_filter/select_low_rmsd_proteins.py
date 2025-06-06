import pandas as pd
import glob
import os

# Set RMSD threshold
rmsd_threshold = 2.0

# Folder with your CSV files
csv_folder = r"C:\Users\Aina Rif'ah\Downloads\UKM\FOR FYPPP\csv file assam result"
csv_files = glob.glob(os.path.join(csv_folder, "*.csv"))

# Storage for selected proteins
selected_rows = []

for file_path in csv_files:
    try:
        # Read CSV and clean column headers
        df = pd.read_csv(file_path, on_bad_lines='skip')
        df.columns = [col.strip() for col in df.columns]
        df["Source File"] = os.path.basename(file_path)

        # Find RMSD column (case-insensitive match)
        rmsd_col = None
        for col in df.columns:
            if "RMSD" in col.upper():
                rmsd_col = col
                break
        if not rmsd_col:
            print(f"⚠️ No RMSD column found in {file_path}")
            continue

        # Extract float from RMSD column safely
        df["RMSD_clean"] = pd.to_numeric(df[rmsd_col].astype(str).str.extract(r"(\d+\.\d+)")[0], errors='coerce')

        # Drop rows with no RMSD
        df = df.dropna(subset=["RMSD_clean"])

        # Filter by low RMSD
        low_rmsd = df[df["RMSD_clean"] <= rmsd_threshold]

        # Take top 2 (lowest RMSD)
        top2 = low_rmsd.nsmallest(2, "RMSD_clean")

        if not top2.empty:
            selected_rows.append(top2)

    except Exception as e:
        print(f"⚠️ Error processing {file_path}: {e}")

# Combine and export
if selected_rows:
    combined_df = pd.concat(selected_rows, ignore_index=True)
    output_file = os.path.join(csv_folder, "Selected_Proteins_LowRMSD.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"\n✅ Done! Saved top 2 proteins per file to:\n{output_file}")
else:
    print("⚠️ No valid proteins found.")
