# @ModuleName: Merge summary tables
# @Function:
# @Author: ggl
# @Time: 2025/10/24 10:24
import pandas as pd
import numpy as np
import subprocess
import os

def clean_filename(name):
    """Remove path and extension, keep only strain name"""
    import os
    base = os.path.basename(name)  # Remove path
    for ext in [".fna", ".fa", ".fasta"]:
        if base.endswith(ext):
            base = base.replace(ext, "")
    return base

# -------------------------------
# User-defined section
# -------------------------------
output_dir = "./output"
padloc_file = os.path.join(output_dir, "padloc_merged_summary.tab")
plasmid_file = os.path.join(output_dir, "all_samples_plasmid_summary.tab")
res_file = os.path.join(output_dir, "all_samples_res_summary.tab")
output_xlsx = os.path.join(output_dir, "merged_summary.xlsx")

padloc_df  = pd.read_csv(padloc_file, sep="\t")
plasmid_df  = pd.read_csv(plasmid_file, sep="\t")
res_df  = pd.read_csv(res_file, sep="\t")

# This is your predefined 139 plasmid replicon headers (example, replace with your full list)
plasmid_headers = [
    "Col(BS512)_1","Col(IMGS31)_1","Col(IRGK)_1","Col(KPHS6)_1","Col(MG828)_1","Col(MP18)_1","Col(VCM04)_1","Col(Ye4449)_1","Col(pHAD28)_1","Col156_1","Col3M_1","Col440II_1","Col440I_1","Col8282_1","ColE10_1","ColKP3_1","ColRNAI_1","ColpEC648","ColpVC_1","IncA_1","IncB/O/K/Z_1","IncB/O/K/Z_2","IncB/O/K/Z_3","IncB/O/K/Z_4","IncC_1","IncFIA(HI1)(pAR0022)_1","IncFIA(HI1)_1","IncFIA(pBK30683)_1","IncFIA_1","IncFIB(AP001918)_1","IncFIB(H89-PhagePlasmid)_1","IncFIB(K)(pCAV1099-114)_1","IncFIB(K)_1","IncFIB(S)_1","IncFIB(pB171)_1","IncFIB(pCTU1)_1","IncFIB(pCTU3)_1","IncFIB(pECLA)_1","IncFIB(pENTAS01)_1","IncFIB(pHCM2)_1","IncFIB(pKPHS1)_1","IncFIB(pLF82-PhagePlasmid)_1","IncFIB(pN55391)_1","IncFIB(pNDM-Mar)_1","IncFIB(pQil)_1","IncFIC(FII)_1","IncFII(29)_1","IncFII(Cf)_1","IncFII(K)_1","IncFII(S)_1","IncFII(SARC14)_1","IncFII(Y)_1","IncFII(Yp)_1","IncFII(p14)_1","IncFII(p96A)_1","IncFII(pAMA1167-NDM-5)_1","IncFII(pAR0022)_1","IncFII(pBK30683)_1","IncFII(pCRY)_1","IncFII(pCTU2)_1","IncFII(pCoo)_1","IncFII(pECLA)_1","IncFII(pEH01)_1","IncFII(pENTA)_1","IncFII(pHN7A8)_1","IncFII(pKP91)_1","IncFII(pKPX1)_1","IncFII(pMET)_1","IncFII(pRSB107)_1","IncFII(pSE11)_1","IncFII(pSFO)_1","IncFII(pYVa12790)_1","IncFII_1","IncHI1A(NDM-CIT)_1","IncHI1A_1","IncHI1B(R27)_1_R27_AF250878","IncHI1B(pNDM-CIT)_1","IncHI1B(pNDM-MAR)_1","IncHI2A_1","IncHI2_1","IncI(Gamma)_1","IncI1-I(Alpha)_1","IncI2(Delta)_1","IncI2_1","IncL_1","IncM1_1","IncM2_1","IncN2_1","IncN3_1","IncN4_1","IncN_1","IncP1_1","IncP1_2","IncP1_3","IncP1_4","IncP6_1","IncQ1_1","IncQ2_1","IncR_1","IncT_1","IncU_1","IncW_1","IncX10_1","IncX1_1","IncX1_2","IncX1_3","IncX1_4","IncX2_1","IncX3(pEC14)_1","IncX3_1","IncX4_1","IncX4_2","IncX5_1","IncX5_2","IncX6_1","IncX8_1","IncX9_1","IncY_1","p0111_1","pADAP_1","pEC4115_1","pENTAS02_1","pESA2_1","pIP31758(p153)_1","pIP31758(p59)_1","pIP32953_1","pKP1433_1","pKPC-CAV1193_1","pKPC-CAV1320_1","pKPC-CAV1321_1","pSL483_1","pSM22_1","pYE854_1","repA(dmsm701b_NDM1)_1","repA(pENTd4a)_1","repA(pKOX)_1_pKOX_CP026273","repB(R1701)_1","repB_KLEB_VIR_AP006726","repE(pEh60-7)_1"
]

# -------------------------------
# 1️⃣ Read padloc summary (main table)
# -------------------------------
padloc_df = pd.read_csv(padloc_file, sep="\t")
padloc_df = padloc_df.rename(columns={"NUM_FOUND": "Defense systems numbers"})
padloc_df["#FILE"] = padloc_df["#FILE"].astype(str)

# Create main DataFrame with padloc #FILE as index
merged_df = pd.DataFrame({"#FILE": padloc_df["#FILE"]})
merged_df["Species"] = 1

# -------------------------------
# 2️⃣ Process plasmid summary
# -------------------------------
plasmid_df = pd.read_csv(plasmid_file, sep="\t")
plasmid_df = plasmid_df.rename(columns={"NUM_FOUND": "Plasmid replicon numbers"})
plasmid_df["#FILE"] = plasmid_df["#FILE"].astype(str)
plasmid_df["#FILE"] = plasmid_df["#FILE"].apply(clean_filename)

# Convert values: numbers (0-100) → 1, "." → 0
plasmid_numeric = plasmid_df.iloc[:, 2:].replace(".", 0)
plasmid_numeric = plasmid_numeric.applymap(lambda x: 1 if str(x).strip() != "0" and str(x).strip() != "." else 0)

# Keep only the 139 predefined headers (add missing columns as 0)
for col in plasmid_headers:
    if col not in plasmid_numeric.columns:
        plasmid_numeric[col] = 0
plasmid_numeric = plasmid_numeric[plasmid_headers]  # Sort by predefined order

# Merge plasmid number column + replicon boolean columns
plasmid_final = pd.concat([plasmid_df[["#FILE", "Plasmid replicon numbers"]], plasmid_numeric], axis=1)

# Left join to merged_df
merged_df = merged_df.merge(plasmid_final, on="#FILE", how="left")
merged_df.fillna(0, inplace=True)

# -------------------------------
# 3️⃣ Process resistance summary
# -------------------------------
res_df = pd.read_csv(res_file, sep="\t")
res_df = res_df.rename(columns={"NUM_FOUND": "ARGs numbers"})
res_df["#FILE"] = res_df["#FILE"].astype(str)
res_df["#FILE"] = res_df["#FILE"].apply(clean_filename)

res_numeric = res_df.iloc[:, 2:].replace(".", 0)
res_numeric = res_numeric.applymap(lambda x: 1 if str(x).strip() != "0" and str(x).strip() != "." else 0)

res_final = pd.concat([res_df[["#FILE", "ARGs numbers"]], res_numeric], axis=1)
merged_df = merged_df.merge(res_final, on="#FILE", how="left")
merged_df.fillna(0, inplace=True)

# -------------------------------
# 4️⃣ Process padloc summary
# -------------------------------
padloc_df["#FILE"] = padloc_df["#FILE"].apply(clean_filename)
padloc_numeric = padloc_df.iloc[:, 2:].replace(".", 0)
padloc_numeric = padloc_numeric.applymap(lambda x: 1 if str(x).strip() != "0" and str(x).strip() != "." else 0)
padloc_final = pd.concat([padloc_df[["#FILE", "Defense systems numbers"]], padloc_numeric], axis=1)
merged_df = merged_df.merge(padloc_final, on="#FILE", how="left")
merged_df.fillna(0, inplace=True)

# -------------------------------
# 5️⃣ Output Excel file
# -------------------------------
merged_df.to_excel(output_xlsx, index=False)
print(f"✅ Merged result file generated: {output_xlsx}")

# ----------------------------
# 🚀 Automatically run "load_randomForest_model.py"
# ----------------------------
current_dir = os.path.dirname(os.path.abspath(__file__))  # Get current script directory
merge_script = os.path.join(current_dir, "load_randomForest_model.py")

if os.path.exists(merge_script):
    print("\n==============================")
    print("✅ Start executing load_randomForest_model.py...")
    print("==============================\n")
    subprocess.run(["python", merge_script], check=True)
    print("\n==============================")
    print("✅ load_randomForest_model.py completed!")
    print("==============================\n")
else:
    print(f"❌ load_randomForest_model.py script not found: {merge_script}")
