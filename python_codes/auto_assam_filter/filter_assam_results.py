import pandas as pd
import glob
import os

# Set RMSD threshold
rmsd_threshold = 2.0

# Folder where your CSV files are located
csv_folder = r"C:\Users\Aina Rif'ah\Downloads\UKM\FOR FYPPP\csv file assam result"
csv_files = glob.glob(os.path.join(csv_folder, "*.csv"))

# Storage for the best match per file
best_per_file = []

for file_path in csv_files:
    try:
        df = pd.read_csv(file_path, on_bad_lines='skip')  # Skip broken lines
        df.columns = [col.strip() for col in df.columns]  # Clean headers
        df["Source File"] = os.path.basename(file_path)

        # Identify RMSD column
        rmsd_col_candidates = [col for col in df.columns if "RMSD" in col]
        if not rmsd_col_candidates:
            continue
        rmsd_col = rmsd_col_candidates[0]

        # Extract numeric RMSD values
        df["RMSD_clean"] = pd.to_numeric(df[rmsd_col].astype(str).str.extract(r"(\d+\.\d+)")[0], errors='coerce')

        # Find column containing "HUMAN"
        match_col = None
        for col in df.columns:
            if df[col].astype(str).str.upper().str.contains("HUMAN").any():
                match_col = col
                break
        if match_col is None:
            continue

        # Apply filtering
        filtered = df[
            (df["RMSD_clean"] <= rmsd_threshold) &
            (df[match_col].astype(str).str.upper().str.contains("HUMAN"))
        ]

        if not filtered.empty:
            # Get the best match (lowest RMSD) from this file
            best_row = filtered.sort_values(by="RMSD_clean").iloc[0]
            best_per_file.append(best_row)

    except Exception as e:
        print(f"⚠️ Error processing {file_path}: {e}")

# Combine and save results
if best_per_file:
    combined_df = pd.DataFrame(best_per_file)
    output_file = os.path.join(csv_folder, "Top1_Human_FromEachFile.csv")
    combined_df.to_csv(output_file, index=False)
    print(f"\n✅ Done! 1 best match per file saved to: {output_file}")
else:
    print("⚠️ No valid matches found in any files.")


