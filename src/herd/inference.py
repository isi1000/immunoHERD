"""
herd/inference.py — sequences -> embedding -> model -> scores.

The glue between the CPU side (validation) and the GPU side (embedding), and
the thing that builds the results table. No science of its own.
"""

from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

from . import esm, models
from .config import THRESHOLD, LABEL_POS, LABEL_NEG, SP2ID, CORE_HOSTS
from .validation import validate_records
from .fasta import parse_fasta, parse_fasta_text


def predict(
    records: List[Tuple[str, str]],
    host: str,
    model: str = "integrated",
    threshold: float = THRESHOLD,
    device: Optional[str] = None,
    batch_seqs: int = 8,
) -> pd.DataFrame:
    """
    records : list of (id, seq).
    host    : target species, by key or alias.
    model   : 'integrated' (recommended, 15 heads) or 'individual' (6 core).

    Returns a DataFrame sorted by descending score, with columns
    id, length, host, model, score, prediction, warnings.
    Invalid sequences come last, with a NaN score.
    """
    if model not in ("integrated", "individual"):
        raise ValueError("model must be 'integrated' or 'individual'.")
    if model == "individual" and host not in CORE_HOSTS:
        raise ValueError(f"'individual' is only available for the 6 core hosts: {CORE_HOSTS}")

    valid, invalid = validate_records(records)

    rows = []
    if valid:
        seqs = [r.seq for r in valid]
        ids = [r.id for r in valid]
        X = esm.embed_sequences(seqs, ids=ids, device=device, batch_seqs=batch_seqs)
        if model == "integrated":
            scores = models.score_integrated(X, host)
        else:
            scores = models.score_individual(X, host)
        for r, sc in zip(valid, scores):
            sc = float(sc)
            rows.append({
                "id": r.id,
                "length": r.orig_length,
                "host": host,
                "model": model,
                "score": round(sc, 4),
                "prediction": LABEL_POS if sc >= threshold else LABEL_NEG,
                "warnings": "; ".join(r.warnings),
            })

    for inv in invalid:
        rows.append({
            "id": inv.id, "length": 0, "host": host, "model": model,
            "score": np.nan, "prediction": "ERROR", "warnings": inv.reason,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("score", ascending=False, na_position="last").reset_index(drop=True)
    return df


def predict_fasta(path: str, host: str, **kwargs) -> pd.DataFrame:
    """Shortcut: predict straight from a FASTA file."""
    ids, seqs = parse_fasta(path)
    return predict(list(zip(ids, seqs)), host, **kwargs)


def predict_text(text: str, host: str, **kwargs) -> pd.DataFrame:
    """Shortcut: predict from pasted text (FASTA, or a bare sequence)."""
    ids, seqs = parse_fasta_text(text)
    return predict(list(zip(ids, seqs)), host, **kwargs)