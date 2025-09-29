#!/usr/bin/env python3
import argparse, json
ap = argparse.ArgumentParser(description='Classify CCHFV contigs into S, M, L using BLAST hits and contig length')
ap.add_argument('--blast', required=True)
ap.add_argument('--contig_qc', required=True)
ap.add_argument('--out_tsv', required=True)
ap.add_argument('--out_json', required=True)
a = ap.parse_args()

# read lengths
len_tbl = {}
with open(a.contig_qc) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 2:
            len_tbl[parts[0]] = int(parts[1])

def guess_segment(L, title):
    t = title.lower()
    if 'segment s' in t or 'small segment' in t: return 'S'
    if 'segment m' in t or 'medium segment' in t: return 'M'
    if 'segment l' in t or 'large segment' in t: return 'L'
    if L < 2500: return 'S'
    if L < 8000: return 'M'
    return 'L'

# parse BLAST table and score hits
seg_hits = {'S': [], 'M': [], 'L': []}
with open(a.blast) as f:
    next(f)
    for line in f:
        parts = line.strip().split('\t')
        contig = parts[0]
        acc = parts[1]
        pident = float(parts[2])
        alen = int(parts[3])
        stitle = parts[-1]
        L = len_tbl.get(contig, 0)
        seg = guess_segment(L, stitle)
        seg_hits[seg].append((contig, acc, pident, alen, stitle))

selected = {}
for seg, hits in seg_hits.items():
    if not hits:
        selected[seg] = None
        continue
    hits_sorted = sorted(hits, key=lambda x: x[2]*x[3], reverse=True)
    selected[seg] = hits_sorted[0]

keyword_map = {'matin': 'Asia-1', 'afg-09-2990': 'Asia-1', 'semunya': 'Africa-2', 'kosova hoti': 'Europe-1'}

rows = []
for seg, hit in selected.items():
    if hit is None:
        rows.append([seg, '', '', '', '', ''])
        continue
    contig, acc, pident, alen, stitle = hit
    clade = ''
    stl = stitle.lower()
    for k, v in keyword_map.items():
        if k in stl:
            clade = v
            break
    rows.append([seg, acc, f'{pident:.2f}', str(alen), stitle, clade])

with open(a.out_tsv, 'w') as out:
    out.write('segment\tacc\tpident\talen\tstitle\tprovisional_clade\n')
    for r in rows:
        out.write('\t'.join(r) + '\n')

with open(a.out_json, 'w') as j:
    json.dump({'rows': rows}, j, indent=2)
print('Wrote', a.out_tsv, 'and', a.out_json)
