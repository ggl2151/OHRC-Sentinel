#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @ModuleName: OHRC-Sentinel cluster analysis
# @Author: ggl

import os
import sys
import glob
import argparse
import subprocess
import tempfile
import shutil

import numpy as np
import pandas as pd


# ============================================================
# Utility functions
# ============================================================

def run_command(cmd, description):

    print("\n--------------------------------")
    print(f"▶ {description}")
    print("--------------------------------")

    print(" ".join(map(str, cmd)))

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.stdout:
        print(result.stdout)

    if result.returncode != 0:

        if result.stderr:
            print(result.stderr)

        raise RuntimeError(
            f"❌ Command failed: {description}"
        )

    return result


def clean_id(x):

    x = os.path.basename(str(x))

    for ext in [
        ".fna",
        ".fa",
        ".fasta",
        ".fna.gz",
        ".fa.gz",
        ".fasta.gz"
    ]:

        if x.endswith(ext):

            x = x[:-len(ext)]

    return x


def pairwise_snp_distance(sequences):

    names = list(sequences.keys())

    rows = []

    for i in range(len(names)):

        for j in range(i + 1, len(names)):

            name1 = names[i]
            name2 = names[j]

            seq1 = sequences[name1]
            seq2 = sequences[name2]

            valid = [
                a != "-"
                and b != "-"
                and a.upper() != "N"
                and b.upper() != "N"
                for a, b in zip(seq1, seq2)
            ]

            if not any(valid):

                distance = np.nan

            else:

                distance = sum(
                    a.upper() != b.upper()
                    for a, b, v
                    in zip(seq1, seq2, valid)
                    if v
                )

            rows.append({
                "Sample_1": name1,
                "Sample_2": name2,
                "SNP_distance": distance
            })

    return pd.DataFrame(rows)


def read_fasta_alignment(path):

    sequences = {}

    name = None
    seq = []

    with open(path, "r") as f:

        for line in f:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    sequences[name] = "".join(seq)

                name = line[1:].split()[0]

                seq = []

            else:

                seq.append(line)

    if name is not None:

        sequences[name] = "".join(seq)

    return sequences


# ============================================================
# FastBAPS
# ============================================================

def run_fastbaps(core_alignment, output_tsv):

    r_code = r'''
args <- commandArgs(trailingOnly=TRUE)

fasta <- args[1]
outfile <- args[2]

library(fastbaps)
library(ape)

sparse.data <- import_fasta_sparse_nt(
    fasta,
    check.fasta = TRUE
)

sparse.data <- optimise_prior(
    sparse.data,
    type = "optimise.symmetric"
)

baps.hc <- fast_baps(
    sparse.data
)

clusters <- best_baps_partition(
    sparse.data,
    baps.hc
)

result <- data.frame(
    Sample_ID = colnames(sparse.data$snp.matrix),
    FastBAPS_cluster = as.character(clusters),
    stringsAsFactors = FALSE
)

write.table(
    result,
    file = outfile,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
'''

    with tempfile.NamedTemporaryFile(
        mode="w",
        suffix=".R",
        delete=False
    ) as rfile:

        rfile.write(r_code)
        r_script = rfile.name

    try:

        run_command(
            [
                "Rscript",
                r_script,
                core_alignment,
                output_tsv
            ],
            "FastBAPS clustering"
        )

    finally:

        if os.path.exists(r_script):

            os.remove(r_script)


# ============================================================
# MOB-suite
# ============================================================

def run_mobtyper(genomes, output_dir):

    mob_dir = os.path.join(
        output_dir,
        "mobtyper"
    )

    os.makedirs(
        mob_dir,
        exist_ok=True
    )

    sample_results = []

    for genome in genomes:

        sample_id = clean_id(
            genome
        )

        sample_out = os.path.join(
            mob_dir,
            f"{sample_id}_mobtyper.txt"
        )

        run_command(
            [
                "mob_typer",
                "--infile",
                genome,
                "--out_file",
                sample_out
            ],
            f"MOB-typer: {sample_id}"
        )

        if not os.path.exists(sample_out):

            print(
                f"⚠️ MOB-typer output not found: "
                f"{sample_out}"
            )

            continue

        try:

            df = pd.read_csv(
                sample_out,
                sep="\t"
            )

        except Exception as e:

            print(
                f"⚠️ Cannot read MOB-typer output "
                f"for {sample_id}: {e}"
            )

            continue

        if df.empty:

            continue

        # One genome can contain multiple plasmids.
        # Aggregate first to one row per genome.
        if "predicted_mobility" in df.columns:

            mobility = (
                df["predicted_mobility"]
                .astype(str)
                .str.strip()
            )

            total = len(mobility)

            conjugative = (
                mobility.str.lower()
                == "conjugative"
            ).sum()

            mobilizable = (
                mobility.str.lower()
                == "mobilizable"
            ).sum()

            non_mobilizable = (
                mobility.str.lower()
                == "non-mobilizable"
            ).sum()

        else:

            total = 0
            conjugative = 0
            mobilizable = 0
            non_mobilizable = 0

        sample_results.append({

            "Sample_ID": sample_id,

            "MOB_typed_plasmids":
                total,

            "Conjugative_plasmids":
                conjugative,

            "Mobilizable_plasmids":
                mobilizable,

            "Non_mobilizable_plasmids":
                non_mobilizable,

            "Conjugative_proportion":
                (
                    conjugative / total
                    if total > 0
                    else np.nan
                )
        })

    if sample_results:

        mob_summary = pd.DataFrame(
            sample_results
        )

    else:

        mob_summary = pd.DataFrame(
            columns=[
                "Sample_ID",
                "MOB_typed_plasmids",
                "Conjugative_plasmids",
                "Mobilizable_plasmids",
                "Non_mobilizable_plasmids",
                "Conjugative_proportion"
            ]
        )

    mob_summary.to_csv(
        os.path.join(
            output_dir,
            "MOB_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    return mob_summary


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="OHRC-Sentinel cluster-level analysis"
    )

    parser.add_argument(
        "--genomes",
        required=True,
        help="Genome directory or glob pattern"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Metadata TSV file"
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="merged_summary_with_prediction.xlsx"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output directory"
    )

    parser.add_argument(
        "--cpus",
        type=int,
        default=8
    )

    parser.add_argument(
        "--min-genomes",
        type=int,
        default=2
    )

    args = parser.parse_args()

    output_dir = os.path.abspath(
        args.output
    )

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    # ========================================================
    # 1. Read metadata
    # ========================================================

    metadata = pd.read_csv(
        args.metadata,
        sep="\t"
    )

    required = {
        "Sample_ID",
        "ST",
        "Source"
    }

    missing = (
        required
        - set(metadata.columns)
    )

    if missing:

        raise ValueError(
            "❌ Metadata missing columns: "
            + ", ".join(sorted(missing))
        )

    metadata["Sample_ID"] = (
        metadata["Sample_ID"]
        .astype(str)
    )

    # ========================================================
    # 2. Find genomes
    # ========================================================

    if os.path.isdir(args.genomes):

        genome_files = sorted(
            glob.glob(
                os.path.join(
                    args.genomes,
                    "*.fna"
                )
            )
        )

    else:

        genome_files = sorted(
            glob.glob(args.genomes)
        )

    if not genome_files:

        raise FileNotFoundError(
            f"❌ No genome files found: "
            f"{args.genomes}"
        )

    genome_map = {
        clean_id(x): x
        for x in genome_files
    }

    metadata = metadata[
        metadata["Sample_ID"]
        .isin(genome_map.keys())
    ].copy()

    if metadata.empty:

        raise ValueError(
            "❌ No metadata Sample_ID matches "
            "the genome filenames."
        )

    # ========================================================
    # 3. Run MOB-suite
    # ========================================================

    mob_summary = run_mobtyper(
        genome_files,
        output_dir
    )

    # ========================================================
    # 4. Load RF prediction summary
    # ========================================================

    summary_df = pd.read_excel(
        args.summary
    )

    summary_df["#FILE"] = (
        summary_df["#FILE"]
        .astype(str)
        .apply(clean_id)
    )

    summary_df = summary_df.rename(
        columns={
            "#FILE": "Sample_ID"
        }
    )

    # ========================================================
    # 5. Merge metadata + summary + MOB
    # ========================================================

    sample_df = metadata.merge(
        summary_df,
        on="Sample_ID",
        how="left"
    )

    sample_df = sample_df.merge(
        mob_summary,
        on="Sample_ID",
        how="left"
    )

    sample_df.to_csv(
        os.path.join(
            output_dir,
            "sample_level_summary.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # 6. Create Snippy working directory
    # ========================================================

    snippy_root = os.path.join(
        output_dir,
        "snippy"
    )

    os.makedirs(
        snippy_root,
        exist_ok=True
    )

    # ========================================================
    # 7. Run Snippy separately for each ST
    # ========================================================

    all_cluster_results = []
    all_pairwise_results = []

    st_groups = (
        metadata
        .groupby("ST")
    )

    for st, group in st_groups:

        group = group.copy()

        if len(group) < args.min_genomes:

            print(
                f"⚠️ {st}: only {len(group)} genomes, "
                "skipping FastBAPS/Snippy."
            )

            continue

        print("\n================================")
        print(
            f"▶ Processing {st}: "
            f"{len(group)} genomes"
        )
        print("================================")

        st_dir = os.path.join(
            snippy_root,
            str(st)
        )

        os.makedirs(
            st_dir,
            exist_ok=True
        )

        # ----------------------------------------------------
        # Select reference
        # ----------------------------------------------------

        reference_id = sorted(
            group["Sample_ID"]
        )[0]

        reference_genome = genome_map[
            reference_id
        ]

        # ----------------------------------------------------
        # Snippy
        # ----------------------------------------------------

        snippy_dirs = []

        for sample_id in sorted(
            group["Sample_ID"]
        ):

            genome = genome_map[
                sample_id
            ]

            outdir = os.path.join(
                st_dir,
                f"snippy_{sample_id}"
            )

            snippy_dirs.append(
                outdir
            )

            if os.path.exists(
                os.path.join(
                    outdir,
                    "snps.vcf"
                )
            ):

                print(
                    f"✓ Snippy already exists: "
                    f"{sample_id}"
                )

                continue

            run_command(
                [
                    "snippy",
                    "--cpus",
                    str(args.cpus),
                    "--outdir",
                    outdir,
                    "--ref",
                    reference_genome,
                    "--ctgs",
                    genome
                ],
                f"Snippy: {st} / {sample_id}"
            )

        # ----------------------------------------------------
        # snippy-core
        # ----------------------------------------------------

        core_prefix = os.path.join(
            st_dir,
            "core"
        )

        core_alignment = (
            core_prefix
            + ".aln"
        )

        if not os.path.exists(
            core_alignment
        ):

            run_command(
                [
                    "snippy-core",
                    "--prefix",
                    core_prefix
                ]
                + snippy_dirs,
                f"snippy-core: {st}"
            )

        if not os.path.exists(
            core_alignment
        ):

            print(
                f"⚠️ core.aln not generated "
                f"for {st}"
            )

            continue

        # ----------------------------------------------------
        # FastBAPS
        # ----------------------------------------------------

        fastbaps_file = os.path.join(
            st_dir,
            "fastbaps_clusters.tsv"
        )

        run_fastbaps(
            core_alignment,
            fastbaps_file
        )

        cluster_df = pd.read_csv(
            fastbaps_file,
            sep="\t"
        )

        cluster_df["ST"] = st

        # ----------------------------------------------------
        # Pairwise SNP distances
        # ----------------------------------------------------

        sequences = read_fasta_alignment(
            core_alignment
        )

        pairwise_df = pairwise_snp_distance(
            sequences
        )

        pairwise_df["ST"] = st

        pairwise_df.to_csv(
            os.path.join(
                st_dir,
                "pairwise_SNP_distances.tsv"
            ),
            sep="\t",
            index=False
        )

        all_pairwise_results.append(
            pairwise_df
        )

        # ----------------------------------------------------
        # Add cluster information
        # ----------------------------------------------------

        all_cluster_results.append(
            cluster_df
        )

    # ========================================================
    # 8. Combine FastBAPS results
    # ========================================================

    if all_cluster_results:

        cluster_assignments = pd.concat(
            all_cluster_results,
            ignore_index=True
        )

    else:

        cluster_assignments = pd.DataFrame(
            columns=[
                "Sample_ID",
                "FastBAPS_cluster",
                "ST"
            ]
        )

    cluster_assignments.to_csv(
        os.path.join(
            output_dir,
            "all_fastbaps_clusters.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # 9. Combine pairwise SNP results
    # ========================================================

    if all_pairwise_results:

        pairwise_all = pd.concat(
            all_pairwise_results,
            ignore_index=True
        )

    else:

        pairwise_all = pd.DataFrame(
            columns=[
                "Sample_1",
                "Sample_2",
                "SNP_distance",
                "ST"
            ]
        )

    pairwise_all.to_csv(
        os.path.join(
            output_dir,
            "pairwise_SNP_distances.tsv"
        ),
        sep="\t",
        index=False
    )

    # ========================================================
    # 10. Add FastBAPS clusters to sample-level data
    # ========================================================

    sample_df = sample_df.merge(
        cluster_assignments[
            [
                "Sample_ID",
                "FastBAPS_cluster"
            ]
        ],
        on="Sample_ID",
        how="left"
    )

    # ========================================================
    # 11. Identify ARG and replicon columns
    # ========================================================

    replicon_columns = [
        x
        for x in summary_df.columns
        if x.startswith(
            (
                "Col",
                "Inc",
                "p",
                "repA",
                "repB",
                "repE"
            )
        )
    ]

    exclude_columns = {
        "Sample_ID",
        "Species",
        "Plasmid replicon numbers",
        "ARGs numbers",
        "Defense systems numbers",
        "Prediction_plasmid_numbers"
    }

    replicon_columns = [
        x
        for x in replicon_columns
        if x not in exclude_columns
    ]

    # ARG columns:
    # all columns after ARGs numbers until non-ARG
    arg_start = None

    columns_all = list(
        summary_df.columns
    )

    if "ARGs numbers" in columns_all:

        arg_start = (
            columns_all.index(
                "ARGs numbers"
            )
            + 1
        )

    if arg_start is not None:

        arg_columns = [
            x
            for x in columns_all[arg_start:]
            if x not in [
                "Defense systems numbers"
            ]
        ]

    else:

        arg_columns = []

    # Remove any accidental non-feature columns
    arg_columns = [
        x
        for x in arg_columns
        if x not in [
            "Prediction_plasmid_numbers"
        ]
    ]

    # ========================================================
    # 12. Cluster-level summary
    # ========================================================

    cluster_records = []

    grouped = sample_df.dropna(
        subset=[
            "ST",
            "FastBAPS_cluster"
        ]
    ).groupby(
        [
            "ST",
            "FastBAPS_cluster"
        ]
    )

    for (
        st,
        cluster
    ), group in grouped:

        record = {

            "ST":
                st,

            "FastBAPS_cluster":
                cluster,

            "Genome_count":
                len(group)
        }

        # ----------------------------------------------------
        # Mean predicted plasmid number
        # ----------------------------------------------------

        if (
            "Prediction_plasmid_numbers"
            in group.columns
        ):

            record[
                "Mean_predicted_plasmid_number"
            ] = pd.to_numeric(
                group[
                    "Prediction_plasmid_numbers"
                ],
                errors="coerce"
            ).mean()

        # ----------------------------------------------------
        # Mean observed plasmid replicon number
        # ----------------------------------------------------

        if (
            "Plasmid replicon numbers"
            in group.columns
        ):

            record[
                "Mean_plasmid_replicon_number"
            ] = pd.to_numeric(
                group[
                    "Plasmid replicon numbers"
                ],
                errors="coerce"
            ).mean()

        # ----------------------------------------------------
        # ARG number
        # ----------------------------------------------------

        if "ARGs numbers" in group.columns:

            record[
                "Mean_ARG_number"
            ] = pd.to_numeric(
                group[
                    "ARGs numbers"
                ],
                errors="coerce"
            ).mean()

        # ----------------------------------------------------
        # Source counts
        # ----------------------------------------------------

        source_counts = (
            group["Source"]
            .value_counts()
        )

        for source, count in (
            source_counts.items()
        ):

            record[
                f"{source}_count"
            ] = int(count)

            record[
                f"{source}_percentage"
            ] = (
                count
                / len(group)
                * 100
            )

        # ----------------------------------------------------
        # MOB-suite
        # ----------------------------------------------------

        if (
            "MOB_typed_plasmids"
            in group.columns
        ):

            total_typed = pd.to_numeric(
                group[
                    "MOB_typed_plasmids"
                ],
                errors="coerce"
            ).sum()

            total_conjugative = pd.to_numeric(
                group[
                    "Conjugative_plasmids"
                ],
                errors="coerce"
            ).sum()

            record[
                "MOB_typed_plasmids"
            ] = total_typed

            record[
                "Conjugative_plasmids"
            ] = total_conjugative

            record[
                "Conjugative_proportion"
            ] = (
                total_conjugative
                / total_typed
                if total_typed > 0
                else np.nan
            )

        # ----------------------------------------------------
        # Mean pairwise SNP distance
        # ----------------------------------------------------

        sample_ids = set(
            group["Sample_ID"]
        )

        pairwise_subset = pairwise_all[
            pairwise_all["Sample_1"].isin(
                sample_ids
            )
            &
            pairwise_all["Sample_2"].isin(
                sample_ids
            )
        ]

        if not pairwise_subset.empty:

            record[
                "Mean_pairwise_SNP_distance"
            ] = pairwise_subset[
                "SNP_distance"
            ].mean()

        else:

            record[
                "Mean_pairwise_SNP_distance"
            ] = np.nan

        cluster_records.append(
            record
        )

    cluster_summary = pd.DataFrame(
        cluster_records
    )

    # ========================================================
    # 13. ARG summary
    # ========================================================

    arg_records = []

    for (
        st,
        cluster
    ), group in grouped:

        for arg in arg_columns:

            values = pd.to_numeric(
                group[arg],
                errors="coerce"
            ).fillna(0)

            n_positive = int(
                (values > 0).sum()
            )

            if n_positive == 0:
                continue

            arg_records.append({

                "ST":
                    st,

                "FastBAPS_cluster":
                    cluster,

                "ARG":
                    arg,

                "Genome_count":
                    len(group),

                "Positive_genomes":
                    n_positive,

                "Prevalence_percent":
                    n_positive
                    / len(group)
                    * 100
            })

    arg_summary = pd.DataFrame(
        arg_records
    )

    # ========================================================
    # 14. Plasmid replicon summary
    # ========================================================

    replicon_records = []

    for (
        st,
        cluster
    ), group in grouped:

        for replicon in replicon_columns:

            values = pd.to_numeric(
                group[replicon],
                errors="coerce"
            ).fillna(0)

            n_positive = int(
                (values > 0).sum()
            )

            if n_positive == 0:
                continue

            replicon_records.append({

                "ST":
                    st,

                "FastBAPS_cluster":
                    cluster,

                "Plasmid_replicon":
                    replicon,

                "Genome_count":
                    len(group),

                "Positive_genomes":
                    n_positive,

                "Prevalence_percent":
                    n_positive
                    / len(group)
                    * 100
            })

    replicon_summary = pd.DataFrame(
        replicon_records
    )

    # ========================================================
    # 15. Save Excel workbook
    # ========================================================

    output_excel = os.path.join(
        output_dir,
        "OHRC_Sentinel_cluster_summary.xlsx"
    )

    with pd.ExcelWriter(
        output_excel,
        engine="openpyxl"
    ) as writer:

        cluster_summary.to_excel(
            writer,
            sheet_name="Cluster_summary",
            index=False
        )

        arg_summary.to_excel(
            writer,
            sheet_name="ARG_summary",
            index=False
        )

        replicon_summary.to_excel(
            writer,
            sheet_name="Replicon_summary",
            index=False
        )

        sample_df.to_excel(
            writer,
            sheet_name="Sample_clusters",
            index=False
        )

    print("\n========================================")
    print("      OHRC-Sentinel analysis finished")
    print("========================================")
    print(
        f"✅ Cluster summary:\n{output_excel}"
    )
    print(
        f"✅ FastBAPS clusters:\n"
        f"{os.path.join(output_dir, 'all_fastbaps_clusters.tsv')}"
    )
    print(
        f"✅ Pairwise SNP distances:\n"
        f"{os.path.join(output_dir, 'pairwise_SNP_distances.tsv')}"
    )
    print(
        f"✅ MOB summary:\n"
        f"{os.path.join(output_dir, 'MOB_summary.tsv')}"
    )
    print(
        f"✅ Sample-level summary:\n"
        f"{os.path.join(output_dir, 'sample_level_summary.tsv')}"
    )


if __name__ == "__main__":

    main()