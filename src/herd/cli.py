"""
herd/cli.py — command line interface.

    herd predict proteins.fasta --host Bos_taurus --out results.csv
    herd hosts

No science here: it turns arguments into calls to `herd.inference`.
"""

import argparse
import sys
from pathlib import Path

from .config import SPECIES_ORDER, CORE_HOSTS, THRESHOLD, LABEL_POS


# --- herd hosts --------------------------------------------------------------
def _cmd_hosts(args) -> int:
    core = set(CORE_HOSTS)
    print("Hosts supported by the integrated model (15 heads):\n")
    for sp in SPECIES_ORDER:
        print(f"  {'core     ' if sp in core else 'auxiliary'} {sp}")
    print(
        "\ncore       also has its own individual model (--model individual)."
        "\nauxiliary  exists only as a head of the integrated model. These nine"
        "\n           heads are not validated to the same level as the core"
        "\n           hosts: interpret their scores with caution."
    )
    return 0


# --- herd predict ------------------------------------------------------------
def _cmd_predict(args) -> int:
    fasta = Path(args.fasta)
    if not fasta.is_file():
        print(f"herd: no such file '{fasta}'", file=sys.stderr)
        return 2

    out = Path(args.out)
    if out.parent and not out.parent.exists():
        print(f"herd: output directory '{out.parent}' does not exist", file=sys.stderr)
        return 2

    # Deferred import: this pulls in torch and TensorFlow, which are slow.
    from .models import resolve_host
    from .inference import predict_fasta

    sp = resolve_host(args.host)
    if sp is None:
        print(f"herd: unknown host '{args.host}'. "
              f"Run 'herd hosts' to see the list.", file=sys.stderr)
        return 2
    if args.model == "individual" and sp not in CORE_HOSTS:
        print(f"herd: no individual model for '{sp}'. "
              f"Only available for: {', '.join(CORE_HOSTS)}.", file=sys.stderr)
        return 2
    if sp != args.host:
        print(f"herd: '{args.host}' interpreted as '{sp}'.", file=sys.stderr)
    if sp not in CORE_HOSTS:
        print(f"herd: warning — '{sp}' is an auxiliary head, not validated to "
              f"the same level as the core hosts.", file=sys.stderr)

    df = predict_fasta(
        str(fasta),
        host=sp,
        model=args.model,
        threshold=args.threshold,
        device=args.device,
        batch_seqs=args.batch_seqs,
    )

    df.to_csv(out, index=False)

    total = len(df)
    errors = int((df["prediction"] == "ERROR").sum()) if total else 0
    positives = int((df["prediction"] == LABEL_POS).sum()) if total else 0
    print(f"herd: {total} sequences | {positives} above threshold "
          f"{args.threshold} | {errors} with errors", file=sys.stderr)
    print(f"herd: results written to {out}", file=sys.stderr)
    return 0


# --- Entry point ------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    from . import __version__

    p = argparse.ArgumentParser(
        prog="herd",
        description="HERD — host-aware prediction of protein immunogenicity.",
    )
    p.add_argument("--version", action="version", version=f"herd {__version__}")
    sub = p.add_subparsers(dest="command", required=True, metavar="COMMAND")

    pr = sub.add_parser("predict", help="Predict from a FASTA file.")
    pr.add_argument("fasta", help="Input FASTA file.")
    pr.add_argument("--host", required=True,
                    help="Target host (see 'herd hosts').")
    pr.add_argument("--out", required=True,
                    help="Output CSV file (required).")
    pr.add_argument("--model", choices=("integrated", "individual"),
                    default="integrated",
                    help="Model to use (default: integrated).")
    pr.add_argument("--threshold", type=float, default=THRESHOLD,
                    help=f"Decision threshold (default: {THRESHOLD}).")
    pr.add_argument("--batch-seqs", type=int, default=8, dest="batch_seqs",
                    help="Sequences per batch when embedding (default: 8).")
    pr.add_argument("--device", default=None,
                    help="Torch device: cpu or cuda (default: automatic).")
    pr.set_defaults(func=_cmd_predict)

    ho = sub.add_parser("hosts", help="List the available hosts.")
    ho.set_defaults(func=_cmd_hosts)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
