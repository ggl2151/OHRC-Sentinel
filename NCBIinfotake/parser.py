# @ModuleName: parser
# @Function: 
# @Author: ggl
# @Time: 2026/7/21 17:04
FIELD_ALIAS = {

    "host": [
        "host",
        "nat_host",
        "host_name",
        "host_organism",
        "host_scientific_name",
        "specific_host",
        "host_common_name",
        "host_species"
    ],

    "isolation_source": [
        "isolation_source",
        "isolate_source",
        "sample_source",
        "source"
    ],

    "country":[
        "country",
        "geo_loc_name",
        "geographic location"
    ],

    "collection_date":[
        "collection_date",
        "collection date"
    ],

    "strain":[
        "strain",
        "isolate"
    ],

    "organism":[
        "organism",
        "organism name"
    ],

    "host_disease":[
        "host_disease",
        "disease",
        "disease_state"
    ],

    "description":[
        "description",
        "title",
        "sample name",
        "sample_title"
    ]
}

def get_field(meta, field):

    aliases = FIELD_ALIAS.get(field, [])

    for key in aliases:

        if key in meta:

            value = meta[key]

            if value:

                return str(value).strip()

    return ""

def normalize_metadata(meta):

    # 保留所有原始字段
    result = meta.copy()

    # 标准化 host
    result["host"] = (
        get_field(meta, "host")
        or get_field(meta, "nat_host")
        or get_field(meta, "host_name")
        or ""
    )

    # 标准化 isolation_source
    result["isolation_source"] = (
        get_field(meta, "isolation_source")
        or get_field(meta, "isolation-source")
        or ""
    )

    # 标准化 country
    result["country"] = (
        get_field(meta, "country")
        or get_field(meta, "geo_loc_name")
        or ""
    )

    # 其它标准字段
    result["collection_date"] = get_field(meta, "collection_date")

    result["strain"] = get_field(meta, "strain")

    result["organism"] = get_field(meta, "organism")

    result["host_disease"] = get_field(meta, "host_disease")

    result["description"] = get_field(meta, "description")

    # 自动推断 host（保留）
    result = infer_host(result)

    return result


def build_search_text(meta):

    text=[]

    for value in meta.values():

        if value:

            text.append(str(value).lower())

    return " ".join(text)


# ==========================================================
# Infer host from metadata
# ==========================================================

# ==========================================================
# Infer host from all metadata
# ==========================================================

def infer_host(meta):
    """
    根据整个 metadata 推断 host。

    优先级：
        1. 如果已有 host，则直接返回
        2. 扫描整个 metadata
        3. 找到宿主后自动补充 host 字段
    """

    # 已有host，不处理
    if meta.get("host"):
        return meta

    # ------------------------------------------------------
    # 将整个metadata拼接成一个字符串
    # ------------------------------------------------------

    text = " ".join(
        str(v)
        for v in meta.values()
        if v is not None and str(v).strip() != ""
    ).lower()

    # ------------------------------------------------------
    # Human
    # ------------------------------------------------------

    human_keywords = [

        "human",
        "homo sapiens",

        "patient",
        "clinical",

        "blood",
        "urine",
        "stool",
        "feces",
        "faeces",

        "wound",
        "pus",
        "sputum",

        "healthy volunteer"

    ]

    for k in human_keywords:

        if k in text:

            meta["host"] = "Human"

            return meta

    # ------------------------------------------------------
    # Animal
    # ------------------------------------------------------

    animal_map = {

        "Chicken":[
            "gallus gallus",
            "gallus",
            "broiler",
            "layer",
            "hen",
            "chicken"
        ],

        "Duck":[
            "anas platyrhynchos",
            "mallard",
            "anas",
            "duck"
        ],

        "Pig":[
            "sus scrofa",
            "porcine",
            "swine",
            "hog",
            "pig"
        ],

        "Cattle":[
            "bos taurus",
            "bovine",
            "cattle",
            "cow",
            "veal",
            "calf"
        ],

        "Sheep":[
            "ovis",
            "ovine",
            "sheep"
        ],

        "Goat":[
            "capra",
            "goat"
        ],

        "Horse":[
            "equus",
            "horse"
        ],

        "Dog":[
            "canis lupus familiaris",
            "canis",
            "dog"
        ],

        "Cat":[
            "felis catus",
            "felis",
            "cat"
        ],

        "Mouse":[
            "mus musculus",
            "mouse"
        ],

        "Rat":[
            "rattus norvegicus",
            "rattus",
            "rat"
        ],

        "Bird":[
            "wild bird",
            "avian",
            "bird"
        ],

        "Fish":[
            "aquaculture",
            "tilapia",
            "salmon",
            "fish"
        ]
    }

    for animal, keywords in animal_map.items():

        for k in keywords:

            if k in text:

                meta["host"] = animal

                return meta

    # ------------------------------------------------------
    # Environment（不是宿主，只做标记）
    # ------------------------------------------------------

    env_keywords = [

        "wastewater",
        "sewage",
        "river",
        "lake",
        "pond",
        "marine",
        "ocean",
        "soil",
        "sediment",
        "farm",
        "manure",
        "hospital surface"

    ]

    for k in env_keywords:

        if k in text:

            meta["host"] = "Environment"

            return meta

    return meta