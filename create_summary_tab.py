#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @ModuleName: OHRC-Sentinel main pipeline
# @Author: ggl

import os
import sys
import glob
import subprocess
import pandas as pd


def process_result(result, combined_file, header_written, basename):
    """
    Check the subprocess result, remove duplicated headers,
    and append the result to the combined file.
    """

    if result.returncode != 0:
        print(f"❌ {basename} encountered an error:")
        print(result.stderr)
        return None, header_written

    lines = result.stdout.strip().splitlines()

    if not lines:
        print(f"⚠️ No output for {basename}, skipping.")
        return None, header_written

    if header_written and lines[0].startswith("#FILE"):
        lines = lines[1:]
    else:
        header_written = True

    if lines:
        with open(combined_file, "a", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    print(f"✅ Results from {basename} appended to {combined_file}")

    return lines, header_written


# ============================================================
# 1. Command-line arguments
# ============================================================

if len(sys.argv) != 4:
    print(
        "\nUsage:\n"
        "python create_summary_tab.py './*.fna' ./output/ metadata.tsv\n"
    )
    sys.exit(1)

input_pattern = sys.argv[1]
output_dir = os.path.abspath(sys.argv[2])
metadata_file = os.path.abspath(sys.argv[3])

os.makedirs(output_dir, exist_ok=True)

if not os.path.exists(metadata_file):
    raise FileNotFoundError(
        f"❌ Metadata file not found: {metadata_file}"
    )

# ============================================================
# 2. Find genome files
# ============================================================

fna_files = sorted(glob.glob(input_pattern))

if not fna_files:
    print(
        f"❌ No .fna files found matching pattern: {input_pattern}"
    )
    sys.exit(1)

print("\n==============================")
print("      OHRC-Sentinel")
print("==============================")
print(f"Genome files : {len(fna_files)}")
print(f"Output dir   : {output_dir}")
print(f"Metadata     : {metadata_file}")
print("==============================\n")


# ============================================================
# 3. Check metadata
# ============================================================

metadata = pd.read_csv(metadata_file, sep="\t")

required_columns = {"Sample_ID", "ST", "Source"}

missing_columns = required_columns - set(metadata.columns)

if missing_columns:
    raise ValueError(
        "❌ Metadata file is missing required columns: "
        + ", ".join(sorted(missing_columns))
    )

metadata["Sample_ID"] = metadata["Sample_ID"].astype(str)

genome_ids = [
    os.path.splitext(os.path.basename(f))[0]
    for f in fna_files
]

missing_metadata = sorted(
    set(genome_ids) - set(metadata["Sample_ID"])
)

if missing_metadata:
    print("⚠️ The following genome files do not have metadata:")
    for x in missing_metadata:
        print(f"   {x}")

    raise ValueError(
        "❌ Every genome must have a corresponding Sample_ID in metadata.tsv."
    )

# Keep only genomes used in this analysis
metadata = metadata[
    metadata["Sample_ID"].isin(genome_ids)
].copy()

metadata.to_csv(
    os.path.join(output_dir, "metadata_used.tsv"),
    sep="\t",
    index=False
)

# ============================================================
# 4. Output files
# ============================================================

pla_combined_tab = os.path.join(
    output_dir, "all_samples_plasmid.tab"
)

res_combined_tab = os.path.join(
    output_dir, "all_samples_res.tab"
)

pla_summary_tab = os.path.join(
    output_dir, "all_samples_plasmid_summary.tab"
)

res_summary_tab = os.path.join(
    output_dir, "all_samples_res_summary.tab"
)

padlocresult_tab = os.path.join(
    output_dir, "padloc_merged.tab"
)

padlocresult_summary_tab = os.path.join(
    output_dir, "padloc_merged_summary.tab"
)


# Remove previous files
for f in [
    pla_combined_tab,
    res_combined_tab,
    pla_summary_tab,
    res_summary_tab,
    padlocresult_tab,
    padlocresult_summary_tab
]:
    if os.path.exists(f):
        os.remove(f)


# ============================================================
# 5. Run ABRicate and PADLOC
# ============================================================

header_written1 = False
header_written2 = False

for fna in fna_files:

    basename = os.path.splitext(
        os.path.basename(fna)
    )[0]

    print("\n--------------------------------")
    print(f"🚀 Processing: {basename}")
    print("--------------------------------")

    # --------------------------------------------------------
    # PlasmidFinder
    # --------------------------------------------------------

    cmd1 = [
        "abricate",
        "--db",
        "plasmidfinder",
        "--mincov",
        "90",
        "--minid",
        "90",
        fna
    ]

    result1 = subprocess.run(
        cmd1,
        capture_output=True,
        text=True
    )

    _, header_written1 = process_result(
        result1,
        pla_combined_tab,
        header_written1,
        basename
    )

    # --------------------------------------------------------
    # ResFinder
    # --------------------------------------------------------

    cmd2 = [
        "abricate",
        "--db",
        "resfinder",
        "--mincov",
        "90",
        "--minid",
        "90",
        fna
    ]

    result2 = subprocess.run(
        cmd2,
        capture_output=True,
        text=True
    )

    _, header_written2 = process_result(
        result2,
        res_combined_tab,
        header_written2,
        basename
    )

    # --------------------------------------------------------
    # PADLOC
    # --------------------------------------------------------

    padloc_outdir = os.path.join(
        output_dir,
        f"{basename}_padloc"
    )

    os.makedirs(
        padloc_outdir,
        exist_ok=True
    )

    cmd3 = [
        "padloc",
        "--fna",
        fna,
        "--cpu",
        "8",
        "-o",
        padloc_outdir
    ]

    result3 = subprocess.run(
        cmd3,
        text=True
    )

    if result3.returncode != 0:
        print(
            f"⚠️ PADLOC failed for {basename}"
        )


# ============================================================
# 6. Generate ABRicate summary
# ============================================================

print("\n📄 Generating ABRicate summary files...")

subprocess.run(
    [
        "abricate",
        "--summary",
        pla_combined_tab
    ],
    stdout=open(pla_summary_tab, "w"),
    check=True
)

subprocess.run(
    [
        "abricate",
        "--summary",
        res_combined_tab
    ],
    stdout=open(res_summary_tab, "w"),
    check=True
)

print(
    f"✅ Plasmid summary: {pla_summary_tab}"
)

print(
    f"✅ ARG summary: {res_summary_tab}"
)


# ============================================================
# 7. Merge PADLOC results
# ============================================================

columns = [
    "#FILE",
    "SEQUENCE",
    "START",
    "END",
    "STRAND",
    "GENE",
    "COVERAGE",
    "COVERAGE_MAP",
    "GAPS",
    "%COVERAGE",
    "%IDENTITY",
    "DATABASE",
    "ACCESSION",
    "PRODUCT",
    "RESISTANCE"
]

merged_df = pd.DataFrame(columns=columns)

csv_files = sorted(
    glob.glob(
        os.path.join(
            output_dir,
            "*_padloc",
            "*.csv"
        )
    )
)

for csv_file in csv_files:

    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(
            f"⚠️ Cannot read PADLOC file: "
            f"{csv_file}\n{e}"
        )
        continue

    if "system" not in df.columns:
        print(
            f"⚠️ {csv_file} does not contain "
            "'system' column, skipping."
        )
        continue

    file_prefix = os.path.basename(
        os.path.dirname(csv_file)
    ).replace("_padloc", "")

    n = len(df)

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
        "RESISTANCE": [""] * n
    })

    merged_df = pd.concat(
        [merged_df, temp_df],
        ignore_index=True
    )

merged_df.to_csv(
    padlocresult_tab,
    sep="\t",
    index=False
)

print(
    f"✅ PADLOC merged file: "
    f"{padlocresult_tab}"
)


# ============================================================
# 8. Generate PADLOC summary
# ============================================================

subprocess.run(
    [
        "abricate",
        "--summary",
        padlocresult_tab
    ],
    stdout=open(
        padlocresult_summary_tab,
        "w"
    ),
    check=True
)

print(
    f"✅ PADLOC summary: "
    f"{padlocresult_summary_tab}"
)


# ============================================================
# 9. Run Integrated_summary_xlsx.py
# ============================================================

current_dir = os.path.dirname(
    os.path.abspath(__file__)
)

merge_script = os.path.join(
    current_dir,
    "Integrated_summary_xlsx.py"
)

if not os.path.exists(merge_script):
    raise FileNotFoundError(
        f"❌ Integrated_summary_xlsx.py not found: "
        f"{merge_script}"
    )

print("\n==============================")
print("▶ Running Integrated_summary_xlsx.py")
print("==============================")

subprocess.run(
    [
        sys.executable,
        merge_script,
        output_dir,
        metadata_file
    ],
    check=True
)

print("\n==============================")
print("✅ OHRC-Sentinel completed")
print("==============================")