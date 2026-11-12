"""Stage 2b: prove/refute g.f.s directly with the transfer-matrix engine.
t_n = e^T Tq^(n+sh) fq for n >= L0 = max(off, -sh) obeys the char-poly recurrence chi (order D) for n >= L0 + D
(indices >= L0).  c = series of N/Q obeys Q (order q) for n > dN.  d = t - c then obeys (Q chi) d = 0 for
n >= N2 = max(L0 + D + q, dN + 1 + D), using indices >= L0.  Hence checking d = 0 on [L0, N2) (plus the listed
terms for off <= n < L0) proves the g.f.; a nonzero d_n refutes it.
usage: stage2b.py OUT.jsonl A...   (env MAXSTATES MAXWORK)
"""
import sys, os, json
sys.setrecursionlimit(10000)
D_ = '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf/'
sys.path.insert(0, D_)
import engine, families, newfam
from gfparse import parse, series, R
from prove2 import parse_entry
from multiprocessing import Pool
SEQ = '/home/user/corpora/oeisdata/seq'
MAXSTATES = int(os.environ.get('MAXSTATES', '60000')); MAXWORK = int(os.environ.get('MAXWORK', '20000000'))
S1 = {}
for l in open(D_ + 'stage1.jsonl'):
    r = json.loads(l)
    if r['s1'] != 'UNPARSED': S1.setdefault(r['A'], []).append(r)

def build(sp):
    preds = sp.get('preds') or [sp['pred']]
    Ts = []
    for p in preds:
        states, succ, acc = engine.make(sp['w'], sp['m'], sp['ext'], p, sp['transpose'], MAXSTATES, MAXWORK)
        Tq, fq, c0, D = engine.quotient(succ, acc)
        Ts.append((Tq, fq, c0, D, len(states)))
    return Ts

def terms(Ts, nmax):
    tot = None
    for Tq, fq, c0, D, _ in Ts:
        a = engine.terms(Tq, fq, c0, nmax)
        tot = a if tot is None else [x + y for x, y in zip(tot, a)]
    return tot

def work(A):
    d = parse_entry(f'{SEQ}/{A[:4]}/{A}.seq'); off = d.get('O', 1); L = d['terms']
    sps = families.specs(d['N']) + newfam.specs(d['N'])
    rows = []
    if not sps:
        return [dict(A=A, kind='gf', status='skipped', detail='no family', formula=r['gf']) for r in S1[A]]
    tried = []
    for sp in sps:
        sh = sp['shift']
        try:
            Ts = build(sp)
        except engine.TooLarge as e:
            tried.append('%s: too large (%s)' % (sp['family'], e)); continue
        Dt = sum(T[3] for T in Ts)
        L0 = max(off, -sh)
        gfs = [(r, parse(r['gf'])) for r in S1[A]]
        need = max([L0 + Dt + len(g.d) + 2 for r, g in gfs] + [len(g.n) + Dt + 2 for r, g in gfs] + [off + len(L) + 1])
        a = terms(Ts, need + sh + 1)
        val = lambda n: a[n + sh]
        mism = [(off + i, val(off + i), v) for i, v in enumerate(L) if off + i >= L0 and val(off + i) != v]
        if mism:
            tried.append('%s: TERM MISMATCH first n=%d computed %d listed %d (%d bad)' % (sp['family'], mism[0][0], mism[0][1], mism[0][2], len(mism)))
            continue
        det = 'family=%s; D=%d (raw %s); w=%d m=%d transpose=%s shift=%d' % (sp['family'], Dt, '+'.join(str(T[4]) for T in Ts), sp['w'], sp['m'], sp['transpose'], sh)
        for r, g in gfs:
            q, dN = len(g.d) - 1, len(g.n) - 1
            N2 = max(L0 + Dt + q, dN + 1 + Dt)
            c = series(g, N2 + 2)
            bad = next((n for n in range(off, N2) if c[n] != (val(n) if n >= L0 else L[n - off])), None)
            row = dict(A=A, kind='gf-direct', formula=r['gf'], idx=r['idx'], D=Dt)
            if bad is None:
                gp = R(g.n, g.d[:-1] + [g.d[-1] + 1]); cp = series(gp, N2 + 2)
                bp = next((n for n in range(off, N2) if cp[n] != (val(n) if n >= L0 else L[n - off])), None)
                row.update(status='PROVED', detail=det + '; checked n=%d..%d; perturbed: %s' % (off, N2 - 1, 'refuted at n=%s' % bp if bp is not None else 'NOT REFUTED'))
            else:
                tv = val(bad) if bad >= L0 else L[bad - off]
                row.update(status='REFUTED at n=%d' % bad, detail=det + '; true a(%d)=%d, g.f. coefficient %s%s' % (bad, tv, c[bad], ' (within listed terms)' if bad < off + len(L) else ''))
                # offset convention check
                for k in (1, -1):
                    gs = R([0] * k + g.n, g.d) if k > 0 else R(g.n, [0] + g.d) if False else None
                    if gs is None: continue
                    cs = series(gs, N2 + 3)
                    if all(cs[n] == (val(n) if n >= L0 else L[n - off]) for n in range(off, N2 + 1)):
                        row['detail'] += '; x*G(x) matches (offset off by one)'
            rows.append(row)
        return rows
    return [dict(A=A, kind='gf', status='skipped' if all('too large' in t for t in tried) else 'TERM MISMATCH', detail=' | '.join(tried), formula=r['gf']) for r in S1[A]]

if __name__ == '__main__':
    out = open(sys.argv[1], 'a')
    As = sys.argv[2:]
    with Pool(int(os.environ.get('NPROC', '4'))) as pool:
        for rows in pool.imap_unordered(work, As, chunksize=1):
            for r in rows:
                out.write(json.dumps(r) + '\n'); out.flush()
                print(r['A'], r['status'], r['detail'][:160], flush=True)
