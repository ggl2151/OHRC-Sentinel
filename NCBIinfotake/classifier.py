# @ModuleName: classifier
# @Function: 
# @Author: ggl
# @Time: 2026/7/21 17:06
# ============================================
# Human
# ============================================

HUMAN_KEYWORDS = [

    "human",

    "homo sapiens",

    "patient",

    "clinical",

    "healthy volunteer",

    "child",

    "adult",

    "infant",

    "neonate"

]

ANIMAL_KEYWORDS = {

    "Chicken":[

        "chicken",

        "broiler",

        "layer",

        "hen",

        "cock",

        "gallus",

        "gallus gallus",

        "gallus gallus domesticus"

    ],

    "Duck":[

        "duck",

        "anas",

        "mallard",

        "anas platyrhynchos"

    ],

    "Pig":[

        "pig",

        "swine",

        "porcine",

        "sus scrofa",

        "hog"

    ],

    "Cattle":[

        "cow",

        "cattle",

        "bovine",

        "bos taurus",

        "calf",

        "veal"

    ],

    "Sheep":[

        "sheep",

        "ovine",

        "ovis"

    ],

    "Goat":[

        "goat",

        "capra"

    ],

    "Horse":[

        "horse",

        "equus"

    ],

    "Dog":[

        "dog",

        "canis",

        "canis lupus familiaris"

    ],

    "Cat":[

        "cat",

        "felis",

        "felis catus"

    ],

    "Mouse":[

        "mouse",

        "mus musculus"

    ],

    "Rat":[

        "rat",

        "rattus",

        "rattus norvegicus"

    ],

    "Bird":[

        "bird",

        "avian",

        "wild bird"

    ],

    "Fish":[

        "fish",

        "tilapia",

        "salmon",

        "aquaculture"

    ]
}

ENVIRONMENT_KEYWORDS = {

    "Wastewater":[

        "wastewater",

        "sewage",

        "waste water",

        "influent",

        "effluent"

    ],

    "River":[

        "river",

        "stream"

    ],

    "Lake":[

        "lake",

        "pond"

    ],

    "Sea":[

        "sea",

        "marine",

        "ocean"

    ],

    "Groundwater":[

        "groundwater"

    ],

    "Soil":[

        "soil"

    ],

    "Sediment":[

        "sediment",

        "mud"

    ],

    "Farm":[

        "farm",

        "manure",

        "compost"

    ],

    "Hospital":[

        "hospital surface",

        "icu surface"

    ],

    "Food":[

        "food",

        "meat",

        "milk",

        "egg",

        "vegetable"

    ]
}

def contains_keyword(text, keywords):

    text = text.lower()

    for k in keywords:

        if k.lower() in text:

            return True

    return False


def classify_source(text):

    text = text.lower()

    ##################################################
    # Human
    ##################################################

    if contains_keyword(text, HUMAN_KEYWORDS):

        return "Human", "Human"

    ##################################################
    # Animal
    ##################################################

    for animal, keywords in ANIMAL_KEYWORDS.items():

        if contains_keyword(text, keywords):

            return "Animal", animal

    ##################################################
    # Environment
    ##################################################

    for env, keywords in ENVIRONMENT_KEYWORDS.items():

        if contains_keyword(text, keywords):

            return "Environment", env

    ##################################################

    return "Unknown", "Unknown"