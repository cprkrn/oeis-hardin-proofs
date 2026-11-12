"""Stage 2a: g.f.s on entries whose recurrence (or polynomial) is already proved.
True sequence t: listed terms for off <= n < start (verified by the provers), then the proved formula.
Let P be the proved recurrence operator (order p, valid for n >= sP, indices used >= off) and Q the g.f.
denominator (order q, Q(0) != 0, numerator degree dN).  d = t - c satisfies the order p+q recurrence QP d = 0
for n >= N2 = max(sP + q, dN + 1 + p); its indices are >= N2-p-q >= off.  So d == 0 on [off, N2) => d == 0.
"""
import sys, json, csv, math
sys.path.insert(0, '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf')
from gfparse import parse, series
from prove2 import parse_entry, parse_rec, parse_poly
from fractions import Fraction
SEQ = '/home/user/corpora/oeisdata/seq'
D = '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf/'
proved = {}
for f in ('proved.csv', 'proved_new.csv'):
    for r in csv.DictReader(open(D + f)):
        proved[r['A']] = (f, r['formula'])
S1 = [json.loads(l) for l in open(D + 'stage1.jsonl')]
out = open(D + 'stage2a.jsonl', 'w')
PERT = 0
for r in S1:
    A = r['A']
    if A not in proved or r['s1'] == 'UNPARSED': continue
    src, ftxt = proved[A]
    d = parse_entry(f'{SEQ}/{A[:4]}/{A}.seq'); off = d.get('O', 1); L = d['terms']
    res = dict(A=A, kind='gf-from-proved-rec(%s)' % src, formula=r['gf'], idx=r['idx'])
    rec = parse_rec([ftxt])
    if rec:
        co, n0, _ = rec; p = max(co)
        start = max(n0 or 0, off + p, p + 2 if src == 'proved.csv' else 0)
        opP = {0: 1, **{j: -c for j, c in co.items()}}   # sum opP[j] t_{n-j} = 0
        sP = start
        def ext(t, n): return sum(c * t[n - j] for j, c in co.items())
    else:
        pp = parse_poly([ftxt])
        if not pp:
            res.update(status='skipped', detail='could not parse proved formula'); out.write(json.dumps(res) + '\n'); continue
        P, deg, n0, _ = pp; p = deg + 1
        start = n0 if n0 is not None else off
        opP = {i: (-1) ** i * math.comb(p, i) for i in range(p + 1)}
        sP = max(start + p, off + p)
        def ext(t, n, P=P): return P(n)
    if start - off > len(L):
        res.update(status='skipped', detail='too few listed terms'); out.write(json.dumps(res) + '\n'); continue
    g = parse(r['gf']); q = len(g.d) - 1; dN = len(g.n) - 1
    if g.d[0] == 0:
        res.update(status='skipped', detail='Q(0)=0'); out.write(json.dumps(res) + '\n'); continue
    N2 = max(sP + q, dN + 1 + p)
    top = max(N2 + 5, off + len(L) + 2)
    t = {n: L[n - off] for n in range(off, start)}
    for n in range(start, top + 1): t[n] = ext(t, n)
    # part 2: listed terms vs true values
    bad_terms = [(n, L[n - off], t[n]) for n in range(start, off + len(L)) if t[n] != L[n - off]]
    c = series(g, top + 1)
    badc = next((n for n in range(off, N2) if c[n] != t[n]), None)
    below = [n for n in range(off) if c[n] != 0]
    if badc is None:
        res.update(status='PROVED', detail=f'proved formula valid n>={start}; P order {p}, Q order {q}, deg N {dN}; checked n={off}..{N2-1}' + (f'; NOTE g.f. has nonzero coefficients below offset at {below}' if below else ''))
    else:
        res.update(status=f'REFUTED at n={badc}', detail=f'true a({badc})={t[badc]}, g.f. coefficient {c[badc]}' + ('' if badc >= off + len(L) else ' (within listed terms)'))
        from gfparse import R as RR
        gx = RR([0] + g.n, g.d); N2x = max(sP + q, dN + 2 + p); cx = series(gx, N2x + 2)
        for n in range(start, N2x + 2): t.setdefault(n, ext(t, n))
        if all(cx[n] == t[n] for n in range(off, N2x)): res['detail'] += f'; x*G(x) PROVED (checked n={off}..{N2x-1}): g.f. is off by a factor x (offset convention)'
    # sanity: perturb the last denominator coefficient on a sample
    if PERT < 200 and badc is None:
        from gfparse import R as RR
        gp = RR(g.n, g.d[:-1] + [g.d[-1] + 1]); cp = series(gp, top + 1)
        bp = next((n for n in range(off, N2) if cp[n] != t[n]), None)
        res['perturb'] = 'refuted at n=%s' % bp if bp is not None else 'NOT REFUTED'; PERT += 1
    if bad_terms: res['bad_terms'] = [list(map(str, x)) for x in bad_terms]
    out.write(json.dumps(res) + '\n')
