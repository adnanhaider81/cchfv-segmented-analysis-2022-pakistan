#!/usr/bin/env python3
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "config.yaml"


def read_fastq_records(path):
    lines = path.read_text().splitlines()
    if len(lines) % 4 != 0:
        raise SystemExit(f"{path} is not valid four-line FASTQ")
    for i in range(0, len(lines), 4):
        header, seq, plus, qual = lines[i:i + 4]
        if not header.startswith("@"):
            raise SystemExit(f"{path} has a FASTQ record without @ header")
        if plus != "+":
            raise SystemExit(f"{path} has a FASTQ record without + separator")
        if len(seq) != len(qual):
            raise SystemExit(f"{path} has sequence and quality lengths that differ")
    return len(lines) // 4


cfg = yaml.safe_load(CONFIG.read_text())
pairs = cfg.get("pairs", [])
if not pairs:
    raise SystemExit("config/config.yaml must contain at least one demo sample")

for item in pairs:
    for key in ("r1", "r2"):
        path = ROOT / item[key]
        if not path.exists():
            raise SystemExit(f"Missing configured demo FASTQ: {path}")
        count = read_fastq_records(path)
        if count == 0:
            raise SystemExit(f"Configured demo FASTQ has no records: {path}")

print(f"Example config OK: {len(pairs)} sample(s)")
