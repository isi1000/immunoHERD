"""
HERD — Host-aware Estimation of antigen immunogenicity and Ranking through
Deep-learning. Shared core behind the web app, the notebook and the CLI.

    from herd import predict_fasta
    df = predict_fasta("proteins.fasta", host="Bos_taurus", model="integrated")
"""

from .inference import predict, predict_fasta, predict_text
from .config import SPECIES_ORDER, CORE_HOSTS, THRESHOLD

__all__ = [
    "predict", "predict_fasta", "predict_text",
    "SPECIES_ORDER", "CORE_HOSTS", "THRESHOLD",
]

__version__ = "0.1.1"