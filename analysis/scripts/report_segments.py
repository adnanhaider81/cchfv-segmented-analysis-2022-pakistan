#!/usr/bin/env python3
import argparse, json
ap = argparse.ArgumentParser(description='Write human friendly report per sample about segments and clades')
ap.add_argument('--segments_json', required=True)
ap.add_argument('--out_txt', required=True)
a = ap.parse_args()
d = json.load(open(a.segments_json))
sel = d.get('rows', [])
clades = [r[5] for r in sel if len(r) >= 6 and r[5]]
mixed = len(set(clades)) > 1
with open(a.out_txt, 'w') as out:
    out.write('Segment selection summary\n')
    for r in sel:
        out.write('\t'.join(r) + '\n')
    out.write('\nInterpretation\n')
    out.write('Segments map to different provisional clades. This suggests possible reassortment.\n' if mixed else 'Segments map to one clade or clade could not be inferred.\n')
print('Wrote', a.out_txt)
