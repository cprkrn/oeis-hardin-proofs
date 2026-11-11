"""Prove Hardin's 'Empirical' recurrences for fixed-width constant-stress tilings:
'Number of (n+1) X (K+1) 0..M arrays with every 2 X 2 subblock having its diagonal sum differing
from its antidiagonal sum by C[, with no adjacent elements equal] (...)'.
State = one row. a(n) = 1^T T^n 1 over the R reachable rows (T maps their span into itself),
so a(n) obeys a linear recurrence of order <= R; if the empirical recurrence Q holds on R
consecutive n it holds for all larger n (Cayley-Hamilton). Smaller n are checked directly."""
import sys, re, os, glob, json, itertools, collections

SEQ, OUT = sys.argv[1], sys.argv[2]
MAXR = int(os.environ.get('MAXR', '3000'))
NAME = re.compile(r'Number of \(n\+1\) X \((\d+)\+1\) 0\.\.(\d+) arrays with every 2 X 2 subblock having its diagonal sum '
                  r'differing from its antidiagonal sum by (\d+)(, with no adjacent elements equal)? \(constant-stress')
REC = re.compile(r'Empirical:\s*a\(n\)\s*=\s*(.+?)(?:\s+for\s+n\s*>\s*(\d+))?\s*\.?\s*$')
TERM = re.compile(r'([+-]?)\s*(\d*)\s*\*?\s*a\(n-(\d+)\)')

def parse(path):
    d = {'S': ''}
    for line in open(path, encoding='utf-8'):
        t, r = line[:2], line[3:].rstrip('\n')
        if t in ('%S', '%T', '%U'): d['S'] += r.split(' ', 1)[1]
        elif t == '%N': d['N'] = r.split(' ', 1)[1]
        elif t == '%O': d['O'] = int(r.split(' ', 1)[1].split(',')[0])
        elif t == '%F': d.setdefault('F', []).append(r.split(' ', 1)[1])
    d['terms'] = [int(x) for x in d['S'].split(',') if x.strip()]
    return d

def parse_rec(F):
    for f in F or []:
        m = REC.search(f)
        if not m or 'a(n-' not in m.group(1): continue
        body = m.group(1).strip()
        if not body.startswith(('+', '-')): body = '+' + body
        co, used = {}, ''
        for t in TERM.finditer(body):
            co[int(t.group(3))] = co.get(int(t.group(3)), 0) + int(t.group(2) or 1) * (-1 if t.group(1) == '-' else 1)
            used += t.group(0)
        if re.sub(r'\s', '', used) != re.sub(r'\s', '', body): return None
        return co, (int(m.group(2)) if m.group(2) else None), f
    return None

def rows_and_moves(w, M, C, noadj, absdiff):
    rows = [r for r in itertools.product(range(M + 1), repeat=w) if not noadj or all(r[i] != r[i + 1] for i in range(w - 1))]
    rid = {r: i for i, r in enumerate(rows)}
    moves = [[] for _ in rows]
    for i, a in enumerate(rows):
        # build next rows b column by column: (a[j]+b[j+1]) - (a[j+1]+b[j]) = +-C
        partial = [(b0,) for b0 in range(M + 1) if not noadj or b0 != a[0]]
        for j in range(w - 1):
            nxt = []
            for p in partial:
                for s in ((C, -C) if absdiff and C else (C,)):
                    v = s + a[j + 1] + p[j] - a[j]   # diag a[j]+b[j+1] minus anti a[j+1]+b[j] = s
                    if 0 <= v <= M and (not noadj or (v != p[j] and v != a[j + 1])):
                        nxt.append(p + (v,))
            partial = nxt
        moves[i] = [rid[b] for b in partial if b in rid]
    return rows, moves

def count(w, M, C, noadj, absdiff, nmax):
    rows, moves = rows_and_moves(w, M, C, noadj, absdiff)
    vec = {i: 1 for i in range(len(rows))}
    reach = set(vec)
    out = {0: len(rows)}
    for n in range(1, nmax + 1):
        nv = collections.defaultdict(int)
        for i, c in vec.items():
            for j in moves[i]: nv[j] += c
        vec = nv; reach |= set(nv)
        out[n] = sum(vec.values())
    return out, len(reach)

results = []
for p in sorted(glob.glob(SEQ + '/*/*.seq')):
    txt = open(p, encoding='utf-8').read()
    if 'constant-stress' not in txt or 'Empirical' not in txt or 'Hardin' not in txt: continue
    d = parse(p); m = NAME.match(d.get('N', ''))
    if not m: continue
    w, M, C, noadj = int(m.group(1)) + 1, int(m.group(2)), int(m.group(3)), bool(m.group(4))
    rec = parse_rec(d.get('F'))
    A = os.path.basename(p)[:-4]
    if rec is None: continue
    co, n0, ftxt = rec; order = max(co); off = d.get('O', 1)
    if (M + 1) ** w > 200000:
        results.append(dict(A=A, status='too large')); continue
    status = None
    for absdiff in (True, False):
        probe, R = count(w, M, C, noadj, absdiff, len(d['terms']) + off)
        if all(probe.get(i + off) == t for i, t in enumerate(d['terms'])):
            break
    else:
        results.append(dict(A=A, status='TERM MISMATCH')); print(A, 'TERM MISMATCH', flush=True); continue
    if R > MAXR:
        results.append(dict(A=A, status='too many states', R=R)); continue
    lo = max(n0 + 1 if n0 is not None else off + order, off + order)
    nmax = lo + R + 1
    a, _ = count(w, M, C, noadj, absdiff, nmax)
    # recurrence derivation: a(n)=1^T T^n 1 holds for n>=0, so Q(E)a obeys order-R recurrence for n>=order
    start_inf = max(lo, order)
    bad = [n for n in range(lo, start_inf + R + 1) if a[n] != sum(c * a[n - j] for j, c in co.items())]
    status = 'PROVED' if not bad else f'REFUTED at n={bad[0]}'
    results.append(dict(A=A, w=w, M=M, C=C, noadj=noadj, absdiff=absdiff, R=R, order=order, valid_from=lo, status=status, formula=ftxt[:160]))
    print(A, status, f'w={w} M={M} C={C} R={R} order={order}', flush=True)
json.dump(results, open(OUT, 'w'), indent=1)
print(collections.Counter(r['status'].split(' at ')[0] for r in results))
