"""
herd/models.py — Keras model loading and scoring.

individual: `model.predict(X)` -> a single sigmoid.
integrated: a sub-model is rebuilt onto the 'species_heads' layer and the
            host's column is taken from it. This deliberately bypasses the
            Lambda layers ('mask' / 'pick'), which do not survive a reload.

Normalization lives INSIDE the .keras file (a Normalization layer adapted at
training time), so the raw (N, 1280) embedding is fed in unscaled.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import tensorflow as tf

from .config import (
    MODELS_ROOT, INDIVIDUAL_DIR, INTEGRATED_FILE,
    SP2ID, CORE_HOSTS, HOST_ALIASES, SPECIES_ORDER,
)

_individual_cache = {}   # host -> keras.Model
_integrated_cache = {}   # 'heads' -> keras.Model (species_heads output)


# --- Host name resolution ---------------------------------------------------
def resolve_host(host: str) -> Optional[str]:
    """Maps a host name onto one of the model's heads, or None."""
    if host is None:
        return None
    h = str(host).strip()
    if h in SP2ID:
        return h
    if h in HOST_ALIASES and HOST_ALIASES[h] in SP2ID:
        return HOST_ALIASES[h]
    genus = h.lower().split("_")[0]          # loose match on genus
    for sp in SP2ID:
        if sp.lower().startswith(genus):
            return sp
    return None


# --- Weight paths -----------------------------------------------------------
def individual_path(host: str) -> Path:
    if host not in CORE_HOSTS:
        raise ValueError(f"No individual model for '{host}'. Core hosts: {CORE_HOSTS}")
    return MODELS_ROOT / INDIVIDUAL_DIR / f"{host}.keras"


def integrated_path() -> Path:
    return MODELS_ROOT / INTEGRATED_FILE


# --- Loading ----------------------------------------------------------------
def load_individual(host: str):
    """Loads and caches the individual model of a core host."""
    if host not in _individual_cache:
        path = individual_path(host)
        _individual_cache[host] = tf.keras.models.load_model(str(path), compile=False)
    return _individual_cache[host]


def load_integrated_heads():
    """
    Loads the integrated model and returns a sub-model whose output is the
    'species_heads' layer: (N, 15) sigmoids. Cached.
    """
    if "heads" not in _integrated_cache:
        path = str(integrated_path())
        try:
            integ = tf.keras.models.load_model(path, compile=False, safe_mode=False)
        except TypeError:  # Keras versions without 'safe_mode'
            integ = tf.keras.models.load_model(path, compile=False)
        heads = tf.keras.Model(integ.inputs[0], integ.get_layer("species_heads").output)
        _integrated_cache["heads"] = heads
    return _integrated_cache["heads"]


# --- Scoring ----------------------------------------------------------------
def score_individual(X: np.ndarray, host: str) -> np.ndarray:
    model = load_individual(host)
    return model.predict(X, verbose=0).ravel()


def score_integrated(X: np.ndarray, host: str) -> np.ndarray:
    sp = resolve_host(host)
    if sp is None:
        raise ValueError(f"Host '{host}' is not in the model. Available: {SPECIES_ORDER}")
    heads = load_integrated_heads()
    H = heads.predict(X, verbose=0)          # (N, 15)
    return H[:, SP2ID[sp]]
