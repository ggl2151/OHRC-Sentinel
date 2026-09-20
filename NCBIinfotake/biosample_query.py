# -*- coding: utf-8 -*-
"""
BioSample Query – 批量提取元数据字段（AI）+ host_type 分类
输出：accession, biosample, biosample_info, host, isolation_source, country, collection_date, host_type
"""
import pandas as pd
import json
from pathlib import Path
from config import USE_AI, OUTPUT_ENCODING
from ncbi_api import get_biosample, get_biosample_metadata
from ai_classifier import batch_classify_with_ai


def process_accession(accession):
    """只获取 BioSample 号和原始元数据（JSON字符串），不进行任何分类"""
    accession = str(accession).strip()
    result = {
        "accession": accession,
        "biosample": None,
        "biosample_info": "",
        "host": "",
        "isolation_source": "",
        "country": "",
        "collection_date": "",
        "host_type": ""
    }

    biosample = get_biosample(accession)
    result["biosample"] = biosample

    if biosample is None:
        return result

    raw_meta = get_biosample_metadata(biosample)  # dict
    print("nihao")
    print(raw_meta)
    # 转为格式化的 JSON 字符串，保留原始信息
    result["biosample_info"] = json.dumps(raw_meta, ensure_ascii=False, indent=2)
    return result


def batch_query(input_file, output_file):
    print("=" * 60)
    print("BioSample Query – AI Metadata Extraction + host_type")
    print("=" * 60)

    # 读取 accession
    df = pd.read_excel(input_file)
    accession_list = (
        df.iloc[:, 0]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    total = len(accession_list)
    if total == 0:
        print("No accession found.")
        return
    print(f"Total accession : {total}\n")

    results = []          # 最终结果列表
    need_ai_list = []     # 待 AI 处理的样本

    success = 0
    failed = 0

    for i, accession in enumerate(accession_list, start=1):
        print(f"[{i}/{total}] {accession}")
        try:
            row = process_accession(accession)
            results.append(row)

            if row.get("biosample"):
                success += 1
                # 如果启用 AI 且有元数据，则加入待处理队列
                if USE_AI and row.get("biosample_info"):
                    try:
                        metadata = json.loads(row["biosample_info"])
                    except:
                        metadata = {}
                    need_ai_list.append({
                        "accession": accession,
                        "metadata": metadata
                    })
            else:
                failed += 1

            print(f"    BioSample : {row.get('biosample')}")
        except Exception as e:
            failed += 1
            print(f"    ERROR : {e}")
            # 添加错误行
            results.append({
                "accession": accession,
                "biosample": None,
                "biosample_info": "",
                "host": "",
                "isolation_source": "",
                "country": "",
                "collection_date": "",
                "host_type": ""
            })
        print("-" * 60)

    # ========== 批量 AI 处理 ==========
    if USE_AI and need_ai_list:
        print(f"\nBatch AI extracting fields for {len(need_ai_list)} samples...")
        try:
            ai_results = batch_classify_with_ai(need_ai_list)
            ai_map = {res["accession"]: res for res in ai_results if "accession" in res}

            for sample in results:
                acc = sample["accession"]
                if acc in ai_map:
                    ai_res = ai_map[acc]
                    sample["host"] = ai_res.get("host", "")
                    sample["isolation_source"] = ai_res.get("isolation_source", "")
                    sample["country"] = ai_res.get("country", "")
                    sample["collection_date"] = ai_res.get("collection_date", "")
                    sample["host_type"] = ai_res.get("host_type", "")
        except Exception as e:
            print(f"Batch AI failed: {e}")

    # ========== 保存 CSV ==========
    out_df = pd.DataFrame(results)
    # 确保列顺序（现在 8 列）
    col_order = ["accession", "biosample", "biosample_info",
                 "host", "isolation_source", "country", "collection_date", "host_type"]
    out_df = out_df[[c for c in col_order if c in out_df.columns]]
    out_df.to_csv(output_file, index=False, encoding=OUTPUT_ENCODING)

    # 打印摘要
    print("\n" + "=" * 60)
    print("Finished!")
    print("=" * 60)
    print(f"Total accession : {total}")
    print(f"Success         : {success}")
    print(f"Failed          : {failed}")
    # 统计 host_type 分布
    if "host_type" in out_df.columns:
        print("\nhost_type distribution:")
        print(out_df["host_type"].value_counts().to_string())
    print()
    print(f"Output file : {Path(output_file).absolute()}")
    print("=" * 60)


if __name__ == "__main__":
    input_file = "accession.xlsx"
    output_file = "output.csv"
    batch_query(input_file, output_file)