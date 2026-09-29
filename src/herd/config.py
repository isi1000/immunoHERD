"""
herd/config.py — fixed constants, no logic. Any hard number or fixed name used
anywhere in the package lives here.
"""

import os
from pathlib import Path

# --- ESM-2 pipeline: identical at training time and at inference ------------
ESM_MODEL   = "esm2_t33_650M_UR50D"   # fair-esm, not the HuggingFace port
REPR_LAYER  = 33
EMBED_DIM   = 1280
SEQ_LENGTH  = 1022
VALID_AA    = "ACDEFGHIKLMNPQRSTVWY"
UNKNOWN_AA  = "X"
THRESHOLD   = 0.5

# --- Integrated model: the order of its 15 heads. Do not reorder ------------
SPECIES_ORDER = [
    "Bos_taurus", "Gallus_gallus", "Homo_sapiens", "Equus_caballus",
    "Sus_scrofa", "Canis_sp", "Camelidae", "Capra_hircus", "Cavia_porcellus",
    "Macaca_sp_", "Mus_musculus", "Non_human_primate", "Oryctolagus_cuniculus",
    "Ovis_aries", "Rattus_sp_",
]
SP2ID = {sp: i for i, sp in enumerate(SPECIES_ORDER)}

# Hosts with an INDIVIDUAL model of their own: the 6 core ones
CORE_HOSTS = [
    "Homo_sapiens", "Bos_taurus", "Sus_scrofa",
    "Gallus_gallus", "Canis_sp", "Equus_caballus",
]

# Host name aliases
HOST_ALIASES = {
    "Canis": "Canis_sp", "Canis_familiaris": "Canis_sp",
    "Canis_lupus_familiaris": "Canis_sp", "Cavia": "Cavia_porcellus",
    "Rattus": "Rattus_sp_", "Macaca": "Macaca_sp_", "Homo": "Homo_sapiens",
}

# --- Where the weights live -------------------------------------------------
# Expected layout:
#   <MODELS_ROOT>/integrated.keras
#   <MODELS_ROOT>/individual/<host>.keras
# Resolved in two steps: HERD_MODELS_DIR if set (web app, tests, alternative
# weights), otherwise the folder shipped next to this module. Resolved from
# __file__ and NOT from the working directory: after a pip install the user's
# cwd has no 'models' folder at all.
_ENV_MODELS = os.environ.get("HERD_MODELS_DIR")
MODELS_ROOT = Path(_ENV_MODELS) if _ENV_MODELS else Path(__file__).resolve().parent / "models"
INDIVIDUAL_DIR  = "individual"
INTEGRATED_FILE = "integrated.keras"

# Output labels
LABEL_POS = "probable immunogen"
LABEL_NEG = "probable non immunogen"