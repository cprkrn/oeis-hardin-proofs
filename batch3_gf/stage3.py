"""Part 2 for already-proved entries: extend verified initial terms with the proved recurrence and compare all listed terms."""
import sys, csv, json
sys.path.insert(0, '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf')
from prove2 import parse_entry, parse_rec
SEQ = '/home/user/corpora/oeisdata/seq'
out = open('stage3.jsonl', 'w'); n = 0
for f in ('proved.csv', 'proved_new.csv'):
    for r in csv.DictReader(open(f)):
        A = r['A']; d = parse_entry(f'{SEQ}/{A[:4]}/{A}.seq'); off = d.get('O', 1); L = d['terms']
        rec = parse_rec([r['formula']])
        if not rec: continue
        co, n0, _ = rec; p = max(co); start = max(n0 or 0, off + p, p + 2 if f == 'proved.csv' else 0)
        t = {k: L[k - off] for k in range(off, min(start, off + len(L)))}
        bad = []
        for k in range(start, off + len(L)):
            t[k] = sum(c * t[k - j] for j, c in co.items())
            if t[k] != L[k - off]: bad.append((k, L[k - off], t[k]))
        n += 1
        out.write(json.dumps(dict(A=A, src=f, start=start, nlisted=len(L), bad=bad)) + '\n')
print(n)
