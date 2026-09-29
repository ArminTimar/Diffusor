"""Index a local PDF library without copying copyrighted sources into the repo.

Usage: python scripts/index_literature.py SOURCE --output tmp/literature
Requires the optional pypdf development dependency.
"""
from pathlib import Path
import argparse
import hashlib
import json

from pypdf import PdfReader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, default=Path("tmp/literature"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    inventory = []
    for path in sorted(args.source.rglob("*.pdf")):
        relative = path.relative_to(args.source)
        target = args.output / relative.with_suffix(".txt")
        target.parent.mkdir(parents=True, exist_ok=True)
        entry = {"file": relative.as_posix(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        try:
            reader = PdfReader(path)
            entry["pages"] = len(reader.pages)
            if not target.exists():
                pages = [f"\n--- PDF PAGE {i + 1} ---\n" + (page.extract_text(extraction_mode="layout") or "")
                         for i, page in enumerate(reader.pages)]
                target.write_text("\n".join(pages), encoding="utf-8")
            entry["extracted_characters"] = target.stat().st_size
        except Exception as exc:
            entry["error"] = str(exc)
        inventory.append(entry)
        print(relative, flush=True)
    (args.output / "inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
