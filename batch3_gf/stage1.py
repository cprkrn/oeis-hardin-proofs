"""Stage 1: parse every Empirical g.f. and compare its expansion with all listed terms."""
import sys, json, re
sys.path.insert(0, '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf')
from gfparse import parse, series, extract
from prove2 import parse_entry
SEQ = '/home/user/corpora/oeisdata/seq'
D = '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf/'
out = open(D + 'stage1.jsonl', 'w')
for A in open(D + 'gf_list.txt').read().split():
    d = parse_entry(f'{SEQ}/{A[:4]}/{A}.seq')
    off = d.get('O', 1)
    gfs = [f for f in d['F'] if 'Empirical g.f.' in f]
    for k, f in enumerate(gfs):
        r = dict(A=A, idx=k, off=off, nterms=len(d['terms']), name=d.get('N', ''), fline=f)
        t = extract(f); r['gf'] = t
        try:
            g = parse(t)
            c = series(g, off + len(d['terms']) + 5)
            r['num'] = [str(x) for x in g.n]; r['den'] = [str(x) for x in g.d]
            # check that c_n == 0 below offset
            low = [n for n in range(off) if c[n] != 0]
            mism = [off + i for i, v in enumerate(d['terms']) if c[off + i] != v]
            if not mism and not low: r['s1'] = 'MATCH'
            elif not mism: r['s1'] = 'MATCH_BUT_NONZERO_BELOW_OFFSET'
            else:
                r['s1'] = 'MISMATCH'; r['first_bad'] = mism[0]; r['nbad'] = len(mism)
                r['oeis'] = d['terms'][mism[0] - off]; r['gfval'] = str(c[mism[0]])
                # shifted convention?
                for sh in (-1, 1, -2, 2):
                    if all(0 <= off + i + sh < len(c) and c[off + i + sh] == v for i, v in enumerate(d['terms'])):
                        r['shift_match'] = sh
        except Exception as e:
            r['s1'] = 'UNPARSED'; r['err'] = str(e)[:80]
        out.write(json.dumps(r) + '\n')
