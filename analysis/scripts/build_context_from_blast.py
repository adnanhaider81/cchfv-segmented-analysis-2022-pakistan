#!/usr/bin/env python3
import argparse, os, time
from Bio import Entrez
ap = argparse.ArgumentParser(description='Fetch top N subject accessions from BLAST for context tree')
ap.add_argument('--blast', required=True)
ap.add_argument('--segment', choices=['S','M','L'], required=True)
ap.add_argument('--hits', type=int, default=25)
ap.add_argument('--out_fasta', required=True)
ap.add_argument('--email', required=False)
ap.add_argument('--api_key', required=False)
a = ap.parse_args()

email = a.email or os.getenv('NCBI_EMAIL')
if not email:
    raise SystemExit('Set --email or env NCBI_EMAIL')
Entrez.email = email
api_key = a.api_key or os.getenv('NCBI_API_KEY')
if api_key:
    Entrez.api_key = api_key

accs = []
with open(a.blast) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        title = parts[-1].lower()
        if a.segment == 'S' and ('segment s' in title or 'small segment' in title):
            accs.append(parts[1])
        elif a.segment == 'M' and ('segment m' in title or 'medium segment' in title):
            accs.append(parts[1])
        elif a.segment == 'L' and ('segment l' in title or 'large segment' in title):
            accs.append(parts[1])
accs = list(dict.fromkeys(accs))[:a.hits]

with open(a.out_fasta, 'w') as out:
    for acc in accs:
        h = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
        out.write(h.read()); h.close()
        time.sleep(0.34)
print('Wrote', a.out_fasta)
