# @ModuleName: ncbi_api
# @Function: 
# @Author: ggl
# @Time: 2026/7/21 17:01
from Bio import Entrez
from xml.etree import ElementTree as ET
import time
import re
import requests



from config import (
    EMAIL,
    NCBI_API_KEY ,
    REQUEST_DELAY,
    MAX_RETRY
)

Entrez.email = EMAIL

if NCBI_API_KEY:
    Entrez.api_key = NCBI_API_KEY

def get_accession_type(accession):

    accession = accession.strip()

    if accession.startswith(("GCA_", "GCF_")):
        return "assembly"

    if accession.startswith(("SAMN", "SAME", "SAMD")):
        return "biosample"

    if accession.startswith(("PRJNA", "PRJEB", "PRJDB")):
        return "bioproject"

    if accession.startswith(("SRR", "ERR", "DRR")):
        return "sra"

    # 剩下基本认为是 nucleotide
    return "nucleotide"


def safe_entrez(method, **kwargs):

    last_error = None

    for i in range(MAX_RETRY):

        try:
            handle = method(**kwargs)
            return handle

        except Exception as e:

            last_error = e

            print(
                f"[Retry {i + 1}/{MAX_RETRY}] "
                f"{type(e).__name__}: {e}"
            )

            if i < MAX_RETRY - 1:
                time.sleep(REQUEST_DELAY)

    print(
        f"[FAILED] Entrez request failed after "
        f"{MAX_RETRY} retries: {last_error}"
    )

    return None

def get_biosample_from_assembly(accession):

    handle = safe_entrez(
        Entrez.esearch,
        db="assembly",
        term=accession
    )

    if handle is None:
        return None

    search = Entrez.read(handle)
    handle.close()

    if len(search["IdList"]) == 0:
        return None

    assembly_id = search["IdList"][0]

    handle = safe_entrez(
        Entrez.esummary,
        db="assembly",
        id=assembly_id,
        report="full"
    )

    if handle is None:
        return None

    summary = Entrez.read(handle)
    handle.close()

    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]

    biosample = doc.get("BioSampleAccn")

    return biosample

def get_biosample_from_nucleotide(accession):

    handle = safe_entrez(
        Entrez.efetch,
        db="nuccore",
        id=accession,
        retmode="xml"
    )

    if handle is None:
        return None

    records = Entrez.read(handle)

    handle.close()

    if len(records) == 0:
        return None

    for x in records[0].get("GBSeq_xrefs", []):

        if x["GBXref_dbname"] == "BioSample":

            return x["GBXref_id"]

    return None


def normalize_metadata_key(key):
    """
    统一 BioSample metadata 字段名称
    """

    if not key:
        return None

    key = str(key).lower().strip()

    # 空格和 - 替换成 _
    key = key.replace("-", "_")
    key = key.replace(" ", "_")

    # 多个连续 _ 合并
    key = re.sub(r"_+", "_", key)

    # 去掉首尾 _
    key = key.strip("_")

    return key


def get_ncbi_biosample_metadata(biosample):
    """
    从 NCBI Entrez BioSample 获取 metadata
    主要用于 SAMN 等 NCBI BioSample accession
    """

    handle = safe_entrez(
        Entrez.efetch,
        db="biosample",
        id=biosample,
        retmode="xml"
    )

    if handle is None:
        print(f"[NCBI FAILED] {biosample}")
        return None

    try:

        xml = handle.read()

    except Exception as e:

        print(f"[NCBI READ ERROR] {biosample}: {e}")
        return None

    finally:

        try:
            handle.close()
        except:
            pass

    if not xml:

        print(f"[NCBI EMPTY RESPONSE] {biosample}")
        return {}

    try:

        root = ET.fromstring(xml)

    except Exception as e:

        print(f"[NCBI XML ERROR] {biosample}: {e}")
        return None

    meta = {}

    for attr in root.iter("Attribute"):

        key = attr.get("attribute_name")

        key = normalize_metadata_key(key)

        if not key:
            continue

        value = attr.text or ""

        meta[key] = value.strip()

    # 如果 NCBI 返回空 BioSampleSet
    if len(meta) == 0:

        biosample_nodes = list(root.iter("BioSample"))

        if len(biosample_nodes) == 0:

            print(f"[NCBI NOT FOUND] {biosample}")

    return meta


def get_ebi_biosample_metadata(biosample):
    """
    从 EBI BioSamples API 获取 metadata
    主要用于 SAMEA accession
    """

    url = f"https://www.ebi.ac.uk/biosamples/samples/{biosample}"

    headers = {
        "Accept": "application/json"
    }

    last_error = None

    for i in range(MAX_RETRY):

        try:

            response = requests.get(
                url,
                headers=headers,
                timeout=30
            )

            # 找不到
            if response.status_code == 404:

                print(f"[EBI NOT FOUND] {biosample}")

                return {}

            response.raise_for_status()

            data = response.json()

            meta = {}

            # EBI BioSamples 的主要 metadata 位于 characteristics
            characteristics = data.get(
                "characteristics",
                {}
            )

            for key, value_info in characteristics.items():

                normalized_key = normalize_metadata_key(key)

                if not normalized_key:
                    continue

                # characteristics 通常是：
                # {
                #   "host": [
                #       {"text": "Bos taurus"}
                #   ]
                # }

                values = []

                if isinstance(value_info, list):

                    for item in value_info:

                        if isinstance(item, dict):

                            value = item.get("text")

                            if value is not None:
                                values.append(
                                    str(value).strip()
                                )

                        elif item is not None:

                            values.append(
                                str(item).strip()
                            )

                elif isinstance(value_info, dict):

                    value = value_info.get("text")

                    if value is not None:
                        values.append(
                            str(value).strip()
                        )

                elif value_info is not None:

                    values.append(
                        str(value_info).strip()
                    )

                # 去掉空值
                values = [
                    v for v in values
                    if v
                ]

                if values:

                    # 多个值使用 ;
                    meta[normalized_key] = "; ".join(values)

            print(
                f"[EBI SUCCESS] {biosample}: "
                f"{len(meta)} metadata fields"
            )

            return meta

        except Exception as e:

            last_error = e

            print(
                f"[EBI Retry {i + 1}/{MAX_RETRY}] "
                f"{biosample}: "
                f"{type(e).__name__}: {e}"
            )

            if i < MAX_RETRY - 1:

                time.sleep(REQUEST_DELAY)

    print(
        f"[EBI FAILED] {biosample} "
        f"after {MAX_RETRY} retries: {last_error}"
    )

    return None


def get_biosample_metadata(biosample):
    """
    根据 BioSample accession 自动选择数据库
    """

    if not biosample:

        return {}

    biosample = str(biosample).strip()

    # -----------------------------
    # NCBI BioSample
    # -----------------------------
    if biosample.startswith("SAMN"):

        print(
            f"[BioSample Source] NCBI: {biosample}"
        )

        return get_ncbi_biosample_metadata(
            biosample
        )

    # -----------------------------
    # EBI / ENA BioSample
    # -----------------------------
    elif biosample.startswith("SAMEA"):

        print(
            f"[BioSample Source] EBI: {biosample}"
        )

        return get_ebi_biosample_metadata(
            biosample
        )

    # -----------------------------
    # DDBJ BioSample
    # 先尝试 EBI
    # -----------------------------
    elif biosample.startswith("SAMD"):

        print(
            f"[BioSample Source] SAMD / trying EBI: "
            f"{biosample}"
        )

        meta = get_ebi_biosample_metadata(
            biosample
        )

        if meta is not None and meta != {}:

            return meta

        print(
            f"[BioSample Source] "
            f"Fallback to NCBI: {biosample}"
        )

        return get_ncbi_biosample_metadata(
            biosample
        )

    # -----------------------------
    # 其他类型
    # 先 NCBI，再 EBI
    # -----------------------------
    else:

        print(
            f"[BioSample Source] Unknown: {biosample}"
        )

        meta = get_ncbi_biosample_metadata(
            biosample
        )

        if meta:

            return meta

        print(
            f"[BioSample Source] "
            f"Trying EBI fallback: {biosample}"
        )

        meta = get_ebi_biosample_metadata(
            biosample
        )

        if meta is not None:

            return meta

        return {}

def get_biosample(accession):
    """
    根据 accession 自动获取 BioSample accession

    支持：
        GCA
        GCF
        CP
        NC
        NZ
        OR
        OX
        SAMN
    """

    accession = accession.strip()

    acc_type = get_accession_type(accession)

    if acc_type == "assembly":

        return get_biosample_from_assembly(accession)

    elif acc_type == "biosample":

        return accession

    elif acc_type == "nucleotide":

        return get_biosample_from_nucleotide(accession)

    else:

        return None