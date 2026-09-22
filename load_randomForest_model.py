#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @ModuleName: Random Forest prediction
# @Author: ggl

import os
import sys
import subprocess

import pandas as pd
import joblib

from sklearn.preprocessing import StandardScaler


# ============================================================
# 1. Command-line arguments
# ============================================================

if len(sys.argv) != 3:

    print(
        "\nUsage:\n"
        "python load_randomForest_model.py OUTPUT_DIR METADATA.tsv\n"
    )

    sys.exit(1)


output_dir = os.path.abspath(sys.argv[1])
metadata_file = os.path.abspath(sys.argv[2])

current_dir = os.path.dirname(
    os.path.abspath(__file__)
)

model_save_dir = os.path.join(
    current_dir,
    "saved_models"
)

input_file_path = os.path.join(
    output_dir,
    "merged_summary.xlsx"
)

output_file_path = os.path.join(
    output_dir,
    "merged_summary_with_prediction.xlsx"
)


# ============================================================
# 2. Check files
# ============================================================

if not os.path.exists(input_file_path):

    raise FileNotFoundError(
        f"❌ Input file not found: {input_file_path}"
    )

if not os.path.exists(metadata_file):

    raise FileNotFoundError(
        f"❌ Metadata file not found: {metadata_file}"
    )


# ============================================================
# 3. Load Random Forest model
# ============================================================

model_path = os.path.join(
    model_save_dir,
    "RandomForest_model.joblib"
)

if not os.path.exists(model_path):

    raise FileNotFoundError(
        f"❌ Model file not found: {model_path}"
    )

model = joblib.load(
    model_path
)

print(
    f"✅ Model loaded: {model_path}"
)


# ============================================================
# 4. Load input data
# ============================================================

df = pd.read_excel(
    input_file_path
)

print(
    f"✅ Input file loaded: "
    f"{input_file_path}"
)

print(
    f"Data shape: {df.shape}"
)


# ============================================================
# 5. Extract exactly the same 141 RF features
# ============================================================

if df.shape[1] < 142:

    raise ValueError(
        "❌ merged_summary.xlsx contains fewer than "
        "142 columns. The Random Forest model expects "
        "141 input features after the first column."
    )

sample_ids = df.iloc[:, 0]

X_new = df.iloc[:, 1:142].copy()

print(
    f"RF feature matrix shape: {X_new.shape}"
)

if X_new.shape[1] != 141:

    raise ValueError(
        f"❌ Expected 141 RF features, "
        f"but found {X_new.shape[1]}."
    )


# ============================================================
# 6. Feature scaling
# ============================================================

scaler_path = os.path.join(
    model_save_dir,
    "scaler.joblib"
)

if os.path.exists(scaler_path):

    scaler = joblib.load(
        scaler_path
    )

    print(
        "✅ Using saved scaler: "
        "scaler.joblib"
    )

    X_scaled = scaler.transform(
        X_new
    )

else:

    print(
        "⚠️ scaler.joblib not found. "
        "Fitting a new scaler."
    )

    scaler = StandardScaler()

    X_scaled = scaler.fit_transform(
        X_new
    )


# ============================================================
# 7. Random Forest prediction
# ============================================================

predictions = model.predict(
    X_scaled
)

print(
    "✅ Random Forest prediction completed"
)


# ============================================================
# 8. Add prediction column
# ============================================================

if "Prediction_plasmid_numbers" in df.columns:

    df.drop(
        columns=["Prediction_plasmid_numbers"],
        inplace=True
    )

# Preserve your original position:
# insert after column C
insert_position = 3

df.insert(
    insert_position,
    "Prediction_plasmid_numbers",
    predictions
)


# ============================================================
# 9. Save prediction result
# ============================================================

df.to_excel(
    output_file_path,
    index=False
)

print(
    f"✅ Prediction file generated:\n"
    f"{output_file_path}"
)


# ============================================================
# 10. Run OHRC-Sentinel cluster analysis
# ============================================================

cluster_script = os.path.join(
    current_dir,
    "run_ohrc_cluster_analysis.py"
)

if not os.path.exists(cluster_script):

    raise FileNotFoundError(
        f"❌ run_ohrc_cluster_analysis.py not found:\n"
        f"{cluster_script}"
    )


print("\n==============================")
print("▶ Starting FastBAPS + Snippy + MOB-suite analysis")
print("==============================")


subprocess.run(
    [
        sys.executable,
        cluster_script,

        "--genomes",
        os.path.dirname(
            os.path.abspath(
                # The genomes are not inside output_dir necessarily.
                # The cluster script will infer them from metadata only
                # if --genome-pattern is provided.
                metadata_file
            )
        ),

        "--metadata",
        metadata_file,

        "--summary",
        output_file_path,

        "--output",
        output_dir
    ],
    check=True
)


print("\n==============================")
print("✅ OHRC-Sentinel cluster analysis completed")
print("==============================")