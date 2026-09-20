# -*- coding: utf-8 -*-
import json
from openai import OpenAI
from config import DEEPSEEK_API_KEY, API_BASE, MODEL

client = OpenAI(api_key=DEEPSEEK_API_KEY, base_url=API_BASE)

def batch_classify_with_ai(sample_list):
    """
    批量提取五字段：host, isolation_source, country, collection_date, host_type
    sample_list: [{"accession": str, "metadata": dict}, ...]
    返回: [{"accession": str, "host": str, "isolation_source": str,
            "country": str, "collection_date": str, "host_type": str}, ...]
    """
    if not sample_list:
        return []

    batch_size = 50   # 可调整
    results = []
    for i in range(0, len(sample_list), batch_size):
        batch = sample_list[i:i+batch_size]
        batch_results = _batch_extract(batch)
        results.extend(batch_results)
    return results

def _batch_extract(batch):
    input_data = []
    for item in batch:
        input_data.append({
            "accession": item["accession"],
            "metadata": item["metadata"]
        })
    message = json.dumps(input_data, ensure_ascii=False, indent=2)

    system_prompt = """
You are an expert in extracting metadata from BioSample records and classifying the sample source.

You will receive a JSON array of samples. Each sample has:
- "accession": unique identifier
- "metadata": all BioSample attributes as key-value pairs.

For each sample, extract or infer the following five fields:

1. host: the biological host (e.g., human, pig, chicken, soil, etc.) – if explicitly provided, use it; otherwise infer from other fields like isolation_source, env_material, etc.
2. isolation_source: the source from which the sample was isolated (e.g., blood, soil, wastewater). Use the "isolation_source" field if present, otherwise infer from env_* fields or description.
3. country: the geographic location (country) – use "geo_loc_name" or "country" field.
4. collection_date: the date of collection if available.
5. host_type: classify the sample into exactly one of the following five categories:
   - Human: if the sample is from human (e.g., host is human, or isolation source is human clinical specimen like blood, urine, sputum).
   - Animal: if the sample is from any animal (vertebrate or invertebrate) – e.g., pig, chicken, fish, crustacean, insect, etc. Even if the exact species is unknown, if it originates from an animal body part, secretion, or surface, classify as Animal.
   - Food: if the sample is from food or food-related products (meat, milk, vegetables, fermented food, etc.).
   - Environment: if the sample is from non-living environmental sources (soil, water, sediment, air, plant surfaces, wastewater, etc.) – if it's not clearly Human, Animal, or Food, it's usually Environment.
   - Unknown: if insufficient information or conflicting.

Rules for host_type:
- Priority: if "host" field exists and is a known animal/human, use that.
- If isolation_source mentions "blood", "urine", "sputum" and there is no host, but implies human clinical, classify as Human.
- If isolation_source or env_material mentions "chicken meat", "pork", "milk" → Food.
- If env_material mentions "crustacean shell", "feces" (animal), "tissue" → Animal (even if host is not explicitly given).
- If sample is from "soil", "water", "sediment", "wastewater" → Environment.
- If sample is from laboratory culture collection or type strain without clear origin → Unknown (or Environment if environmental source is clear).

Return a JSON array with the same order as input. Each element must be:
{
    "accession": "...",
    "host": "...",
    "isolation_source": "...",
    "country": "...",
    "collection_date": "...",
    "host_type": "Human|Animal|Food|Environment|Unknown"
}

Return ONLY valid JSON, no markdown.
"""

    response = client.chat.completions.create(
        model=MODEL,
        temperature=0,
        response_format={"type": "json_object"},
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": message}
        ]
    )

    result_text = response.choices[0].message.content
    parsed = json.loads(result_text)

    # 容错处理
    if isinstance(parsed, dict):
        if "results" in parsed:
            parsed = parsed["results"]
        else:
            parsed = [parsed]
    return parsed