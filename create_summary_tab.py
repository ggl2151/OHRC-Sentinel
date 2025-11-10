# @ModuleName: Generate summary tables by running commands
# @Function:
# @Author: ggl
# @Time: 2025/10/24 11:25
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import subprocess
import pandas as pd

def process_result(result, combined_file, header_written, basename):
    """
    Check the subprocess result, extract lines, remove the header
    (only keep header from the first sample), and append to the combined file.
    Returns the updated header_written status.
    """
    if result.returncode != 0:
        print(f"❌ {basename} encountered an error:\n{result.stderr}\n")
        return None, header_written  # Return None to skip

    lines = result.stdout.strip().splitlines()
    if not lines:
        print(f"⚠️ No output for {basename}, skipping.\n")
        return None, header_written

    # Remove header; keep only the first sample’s header
    if header_written and lines[0].startswith("#FILE"):
        lines = lines[1:]
    else:
        header_written = True

    with open(combined_file, "a", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(f"✅ Results from {basename} appended to {combined_file}\n")
    return lines, header_written

# Usage: python run_abricate.py "./*.fna" "./output/"
if len(sys.argv) != 3:
    print("Example: python create_summary_tab.py './*.fna' ./output/")
    sys.exit(1)

input_pattern = sys.argv[1]
output_dir = sys.argv[2]

# Match all .fna files
fna_files = sorted(glob.glob(input_pattern))
if not fna_files:
    print(f"❌ No .fna files found matching pattern: {input_pattern}")
    sys.exit(1)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

pla_combined_tab = os.path.join(output_dir, "all_samples_plasmid.tab")
res_combined_tab = os.path.join(output_dir, "all_samples_res.tab")
pla_summary_tab = os.path.join(output_dir, "all_samples_plasmid_summary.tab")
res_summary_tab = os.path.join(output_dir, "all_samples_res_summary.tab")
padloc_combined_csv = os.path.join(output_dir, "all_samples_padloc.csv")

# Remove old files if they exist
for f in [pla_combined_tab, res_combined_tab, pla_summary_tab, res_summary_tab, padloc_combined_csv]:
    if os.path.exists(f):
        os.remove(f)

padloc_csv_files = []
header_written1 = False
header_written2 = False

for fna in fna_files:
    basename = os.path.splitext(os.path.basename(fna))[0]
    print(f"🚀 Processing: {basename}")

    # Keep original command format
    cmd1 = f"abricate --db plasmidfinder --mincov 90 --minid 90 {fna}"
    cmd2 = f"abricate --db resfinder --mincov 90 --minid 90 {fna}"

    # Capture command outputs
    result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    result2 = subprocess.run(cmd2, shell=True, capture_output=True, text=True)

    # Main loop
    lines1, header_written1 = process_result(result1, pla_combined_tab, header_written1, basename)
    if lines1 is None:
        continue  # Skip this iteration

    lines2, header_written2 = process_result(result2, res_combined_tab, header_written2, basename)
    if lines2 is None:
        continue  # Skip this iteration

    # Run padloc
    padloc_outdir = os.path.join(output_dir, f"{basename}_padloc")
    os.makedirs(padloc_outdir, exist_ok=True)

    # Construct padloc command
    cmd3 = f"padloc --fna {fna} --cpu 8 -o {padloc_outdir}"
    subprocess.run(cmd3, shell=True, text=True)



# Step 2: Generate summary tables
print("📄 Generating summary files ...")
cmd1_2 = f"abricate --summary {pla_combined_tab} > {pla_summary_tab}"
cmd2_2 = f"abricate --summary {res_combined_tab} > {res_summary_tab}"
subprocess.run(cmd1_2, shell=True, check=True)
subprocess.run(cmd2_2, shell=True, check=True)

print(f"🎉 All files processed! Summary file located at: {pla_summary_tab}")
print(f"🎉 All files processed! Summary file located at: {res_summary_tab}")


# ---------------------------
# Merge all padloc CSV files
# ---------------------------
cmd3_2 = f"mv ./output/*/*.csv ./output/"
subprocess.run(cmd3_2, shell=True, check=True)

# Output files
padlocresult_tab = "./output/padloc_merged.tab"
padlocresult_summary_tab = "./output/padloc_merged_summary.tab"

# Fixed header columns
columns = [
    "#FILE", "SEQUENCE", "START", "END", "STRAND", "GENE",
    "COVERAGE", "COVERAGE_MAP", "GAPS", "%COVERAGE", "%IDENTITY",
    "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"
]

# Create empty DataFrame
merged_df = pd.DataFrame(columns=columns)

# Traverse all *_padloc.csv files under output/
csv_files = sorted(glob.glob("./output/*_padloc.csv"))

for csv_file in csv_files:
    # Read CSV file
    df = pd.read_csv(csv_file)

    # Extract "system" column as GENE
    if "system" not in df.columns:
        print(f"⚠️ {csv_file} does not contain 'system' column, skipping.")
        continue

    n = len(df)
    file_prefix = os.path.basename(csv_file).split("_padloc")[0]

    temp_df = pd.DataFrame({
        "#FILE": [file_prefix] * n,
        "SEQUENCE": [""] * n,
        "START": [""] * n,
        "END": [""] * n,
        "STRAND": [""] * n,
        "GENE": df["system"],
        "COVERAGE": [""] * n,
        "COVERAGE_MAP": [""] * n,
        "GAPS": [""] * n,
        "%COVERAGE": [100] * n,
        "%IDENTITY": [100] * n,
        "DATABASE": [""] * n,
        "ACCESSION": [""] * n,
        "PRODUCT": [""] * n,
        "RESISTANCE": [""] * n,
    })

    merged_df = pd.concat([merged_df, temp_df], ignore_index=True)

# Save as tab-delimited file
merged_df.to_csv(padlocresult_tab, sep="\t", index=False)
print(f"✅ Merged file generated: {padlocresult_tab}")

# Generate padloc summary
cmd3_3 = f"abricate --summary {padlocresult_tab} > {padlocresult_summary_tab}"
subprocess.run(cmd3_3, shell=True, check=True)
print(f"🎉 All files processed! Summary file located at: {padlocresult_summary_tab}")




# ----------------------------
# 🚀 Automatically run “Integrated_summary_xlsx.py”
# ----------------------------
current_dir = os.path.dirname(os.path.abspath(__file__))  # Get current script directory
merge_script = os.path.join(current_dir, "Integrated_summary_xlsx.py")

if os.path.exists(merge_script):
    print("\n==============================")
    print("✅ Executing integrated summary merge script...")
    print("==============================\n")
    subprocess.run(["python", merge_script], check=True)
    print("\n==============================")
    print("✅ Integrated summary merge completed!")
    print("==============================\n")
else:
    print(f"❌ Merge script not found: {merge_script}")
