"""Prove (or refute) R. H. Hardin's 'Empirical' linear recurrences for fixed-width
'n X k 0..m arrays with every element equal/unequal to S neighbours' OEIS entries.

Method: rows are generated top to bottom; state = (previous row, current row).
a(n) = number of valid n-row arrays = x^T T^(n-2) y for n >= 2, so a(n) satisfies the
characteristic polynomial of T (order D). If the empirical recurrence q annihilates
a(n) for D consecutive valid n, then q(E)a (which also satisfies the order-D recurrence)
is identically zero -> the empirical recurrence is proved for all n beyond its offset.
"""
import sys, re, os, glob, json, itertools
from fractions import Fraction

SEQ = sys.argv[1]
OUT = sys.argv[2]
MAXSTATES = int(os.environ.get('MAXSTATES', '1100'))

NAME = re.compile(r'Number of n ?X ?(\d+) 0\.\.(\d+) arrays with every element (equal|unequal) to ([\d, or]+?) '
                  r'(king-move|horizontally, vertically or antidiagonally) adjacent elements, with upper left element zero\.$')

def parse_entry(path):
    d = {'S': ''}
    for line in open(path, encoding='utf-8'):
        tag, rest = line[:2], line[3:].rstrip('\n')
        if tag in ('%S', '%T', '%U'):
            d['S'] += rest.split(' ', 1)[1]
        elif tag == '%N':
            d['N'] = rest.split(' ', 1)[1]
        elif tag == '%O':
            d['O'] = int(rest.split(' ', 1)[1].split(',')[0])
        elif tag == '%F':
            d.setdefault('F', []).append(rest.split(' ', 1)[1])
    d['terms'] = [int(x) for x in d['S'].split(',') if x.strip()]
    return d

REC = re.compile(r'Empirical:\s*a\(n\)\s*=\s*(.+?)(?:\s+for\s+n\s*>\s*(\d+))?\s*\.?\s*$')
TERM = re.compile(r'([+-]?)\s*(\d*)\s*\*?\s*a\(n-(\d+)\)')

def parse_rec(F):
    for f in F or []:
        m = REC.search(f)
        if m and 'a(n-' in m.group(1):
            body = m.group(1).strip()
            if not body.startswith(('+', '-')):
                body = '+' + body
            coefs = {}
            consumed = ''
            for t in TERM.finditer(body):
                c = int(t.group(2) or 1) * (-1 if t.group(1) == '-' else 1)
                coefs[int(t.group(3))] = coefs.get(int(t.group(3)), 0) + c
                consumed += t.group(0)
            if re.sub(r'\s', '', consumed) != re.sub(r'\s', '', body):
                return None  # has constant / polynomial terms -> skip
            return coefs, int(m.group(2)) if m.group(2) else None, f
    return None

def neighbours(kind):
    if kind == 'king-move':
        return [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]
    return [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, 1), (1, -1)]  # horiz, vert, antidiagonal (NE, SW)

def build(k, m, eq, allowed, kind):
    nb = neighbours(kind)
    rows = list(itertools.product(range(m + 1), repeat=k))
    def ok_row(above, cur, below):
        for j in range(k):
            cnt = 0
            for di, dj in nb:
                r = above if di == -1 else below if di == 1 else cur
                if r is None:
                    continue
                jj = j + dj
                if 0 <= jj < k:
                    v = r[jj]
                    cnt += (v == cur[j]) if eq else (v != cur[j])
            if cnt not in allowed:
                return False
        return True
    return rows, ok_row

def count_terms(k, m, eq, allowed, kind, nmax):
    rows, ok = build(k, m, eq, allowed, kind)
    # n = 1
    res = {1: sum(1 for r in rows if r[0] == 0 and ok(None, r, None))}
    # state vector over (prev, cur) where prev already validated
    vec = {}
    for a in rows:
        if a[0] != 0:
            continue
        for b in rows:
            if ok(None, a, b):
                vec[(a, b)] = vec.get((a, b), 0) + 1
    if len(vec) == 0:
        return res, 0
    states = set()
    for n in range(2, nmax + 1):
        res[n] = sum(c for (a, b), c in vec.items() if ok(a, b, None))
        nv = {}
        for (a, b), c in vec.items():
            for cnext in rows:
                if ok(a, b, cnext):
                    nv[(b, cnext)] = nv.get((b, cnext), 0) + c
        vec = nv
        states |= set(vec)
        if len(states) > MAXSTATES * 4:
            pass
    return res, len(rows) ** 2

results = []
files = glob.glob(SEQ + '/*/*.seq')
for p in sorted(files):
    txt = open(p, encoding='utf-8').read()
    if 'Hardin' not in txt or 'Empirical' not in txt:
        continue
    d = parse_entry(p)
    m = NAME.match(d.get('N', ''))
    if not m:
        continue
    k, mm = int(m.group(1)), int(m.group(2))
    eq = m.group(3) == 'equal'
    allowed = set(int(x) for x in re.findall(r'\d+', m.group(4)))
    kind = m.group(5)
    rec = parse_rec(d.get('F'))
    if rec is None:
        continue
    coefs, n0, ftxt = rec
    order = max(coefs)
    D = (mm + 1) ** (2 * k)
    if D > MAXSTATES:
        results.append(dict(A=os.path.basename(p)[:-4], k=k, m=mm, status='too many states', D=D))
        continue
    off = d.get('O', 1)
    start = max(n0 + 1 if n0 is not None else order + off, order + 2)
    nmax = start + D + 2
    terms, _ = count_terms(k, mm, eq, allowed, kind, nmax)
    # 1) compare with OEIS terms
    oeis = d['terms']
    mism = [i + off for i, t in enumerate(oeis) if terms.get(i + off) is not None and terms[i + off] != t]
    if mism:
        results.append(dict(A=os.path.basename(p)[:-4], k=k, m=mm, status='TERM MISMATCH (interpretation or erratum)',
                            first_bad=mism[0], D=D))
        continue
    # 2) check recurrence on D consecutive n (a(n) = x^T T^(n-2) y valid for n >= 2; need n-order >= 2)
    lo = max(start, order + 2)
    bad = [n for n in range(max(min(start, lo), order + 1), lo + D) if terms[n] != sum(c * terms[n - j] for j, c in coefs.items())]
    status = 'PROVED' if not bad else f'REFUTED at n={bad[0]}'
    results.append(dict(A=os.path.basename(p)[:-4], k=k, m=mm, status=status, D=D, order=order, valid_from=lo,
                        formula=ftxt[:160]))
    print(results[-1]['A'], status, f'k={k} m={mm} D={D} order={order}', flush=True)

json.dump(results, open(OUT, 'w'), indent=1)
from collections import Counter
print(Counter(r['status'].split(' at ')[0] for r in results))
