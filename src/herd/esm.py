"""
herd/esm.py — ESM-2 loading and embedding. This is the GPU part.

Four things here must not change:
  - the model is loaded ONCE and reused (lazy singleton)
  - .eval() + torch.no_grad()
  - numpy is returned, never CUDA tensors (this matters on ZeroGPU)
  - the pooling recipe `reps[k, 1:tr+1].mean(0)`

Expects sequences that are ALREADY CLEAN (see validation.py).
"""

from typing import List, Optional

import numpy as np
import torch
from esm import pretrained

from .config import ESM_MODEL, REPR_LAYER, SEQ_LENGTH, EMBED_DIM

_esm = None  # module-level cache: (model, alphabet, batch_converter, device)


def get_esm(device: Optional[str] = None):
    """Loads ESM-2 on first call, reuses it afterwards."""
    global _esm
    if _esm is None:
        model, alphabet = pretrained.load_model_and_alphabet(ESM_MODEL)
        model.eval()
        device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        model = model.to(device)
        batch_converter = alphabet.get_batch_converter(SEQ_LENGTH)
        _esm = (model, alphabet, batch_converter, device)
    return _esm


@torch.no_grad()
def embed_sequences(
    seqs: List[str],
    ids: Optional[List[str]] = None,
    device: Optional[str] = None,
    batch_seqs: int = 8,
) -> np.ndarray:
    """
    Clean sequences -> (N, 1280) float32 matrix, in the order they came in.

    `batch_seqs` is the maximum number of sequences per batch.

    Sequences are batched SORTED BY LENGTH and the result is put back in the
    original order. Every batch is padded to its longest sequence, so grouping
    similar lengths removes almost all the padding: on a real proteome (4,570
    sequences, mean 386 residues, 5% at the 1022 cap) this halves the GPU work.

    It does not change the numbers: each mean is taken over that sequence's
    real residues only, and padding never enters it. Nor does it raise peak
    memory: the worst batch is still `batch_seqs` sequences at full length.
    """
    model, alphabet, batch_converter, dev = get_esm(device)
    n = len(seqs)
    X = np.zeros((n, EMBED_DIM), dtype="float32")
    if n == 0:
        return X

    if ids is None:
        ids = [str(i) for i in range(n)]

    # Longest first, on purpose: if this is going to run out of memory it
    # does so on the first batch rather than halfway through.
    orden = sorted(range(n), key=lambda i: len(seqs[i]), reverse=True)

    for st in range(0, n, batch_seqs):
        idx = orden[st:st + batch_seqs]
        chunk = [seqs[i] for i in idx]
        labels = [ids[i] for i in idx]
        _, _, toks = batch_converter(list(zip(labels, chunk)))
        toks = toks.to(dev)
        reps = model(toks, repr_layers=[REPR_LAYER])["representations"][REPR_LAYER].to("cpu")
        for k, (i, s) in enumerate(zip(idx, chunk)):
            tr = min(SEQ_LENGTH, len(s))
            # Exact recipe: drop BOS, exclude EOS and padding.
            # Written to the ORIGINAL position i, not the batch position.
            X[i] = reps[k, 1:tr + 1].mean(0).numpy().astype("float32")
    return X
