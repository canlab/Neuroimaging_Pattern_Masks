#!/usr/bin/env python3
"""Generate OpenAI embeddings for topic-term documents."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path


DEFAULT_MODEL = "text-embedding-3-large"
API_URL = "https://api.openai.com/v1/embeddings"
DOC_PATTERN = re.compile(r"^(\d+)_(.+)\.csv$")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Generate embeddings for Neurosynth topic-term CSV documents."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=repo_root / "topic_terms_csv",
        help="Directory containing numbered topic CSV files.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=repo_root / f"topic_embeddings_{DEFAULT_MODEL}.csv",
        help="Path for the output embedding matrix CSV.",
    )
    parser.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help="OpenAI embedding model name.",
    )
    parser.add_argument(
        "--pause-seconds",
        type=float,
        default=0.0,
        help="Optional pause between API requests.",
    )
    return parser.parse_args()


def load_documents(input_dir: Path) -> list[tuple[int, str, str]]:
    documents: list[tuple[int, str, str]] = []
    for path in input_dir.glob("*.csv"):
        match = DOC_PATTERN.match(path.name)
        if not match:
            continue
        doc_number = int(match.group(1))
        doc_name = path.stem
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            terms = [row[0].strip() for row in reader if row and row[0].strip()]
        if terms and terms[0].lower() == "term":
            terms = terms[1:]
        text = "\n".join(terms)
        documents.append((doc_number, doc_name, text))
    documents.sort(key=lambda item: item[0])
    if not documents:
        raise FileNotFoundError(f"No topic CSV files found in {input_dir}")
    return documents


def fetch_embedding(text: str, model: str, api_key: str) -> list[float]:
    payload = json.dumps({"input": text, "model": model}).encode("utf-8")
    request = urllib.request.Request(
        API_URL,
        data=payload,
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            body = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"OpenAI API request failed ({exc.code}): {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"Network error while calling OpenAI API: {exc}") from exc

    data = body.get("data")
    if not data or "embedding" not in data[0]:
        raise RuntimeError(f"Unexpected OpenAI API response: {body}")
    return data[0]["embedding"]


def write_embeddings_csv(
    output_csv: Path, document_names: list[str], embeddings: list[list[float]]
) -> None:
    embedding_length = len(embeddings[0])
    for embedding in embeddings[1:]:
        if len(embedding) != embedding_length:
            raise ValueError("Embedding vectors do not all have the same length")

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(document_names)
        for row_idx in range(embedding_length):
            writer.writerow([embedding[row_idx] for embedding in embeddings])


def main() -> int:
    args = parse_args()
    api_key = os.environ.get("OPENAI_API_KEY")
    if not api_key:
        print("OPENAI_API_KEY is not set.", file=sys.stderr)
        return 2

    documents = load_documents(args.input_dir)
    embeddings: list[list[float]] = []

    for index, (_, doc_name, text) in enumerate(documents, start=1):
        print(f"[{index}/{len(documents)}] Embedding {doc_name}...", flush=True)
        embeddings.append(fetch_embedding(text, args.model, api_key))
        if args.pause_seconds > 0:
            time.sleep(args.pause_seconds)

    write_embeddings_csv(
        args.output_csv,
        [doc_name for _, doc_name, _ in documents],
        embeddings,
    )
    print(f"Wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
