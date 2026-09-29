"""
herd/fasta.py — reading sequences. parse_fasta takes a file, parse_fasta_text
takes a pasted string.
"""

from typing import List, Tuple


def parse_fasta(path: str) -> Tuple[List[str], List[str]]:
    """Reads a FASTA file -> (ids, seqs)."""
    ids, seqs, name, chunks = [], [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    ids.append(name)
                    seqs.append("".join(chunks))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        ids.append(name)
        seqs.append("".join(chunks))
    return ids, seqs


def parse_fasta_text(text: str) -> Tuple[List[str], List[str]]:
    """Same, from a string. With no '>', the whole text is one sequence."""
    text = text.strip()
    if not text:
        return [], []
    if not text.lstrip().startswith(">"):
        return ["sequence_1"], ["".join(text.split())]

    ids, seqs, name, chunks = [], [], None, []
    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                ids.append(name)
                seqs.append("".join(chunks))
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if name is not None:
        ids.append(name)
        seqs.append("".join(chunks))
    return ids, seqs