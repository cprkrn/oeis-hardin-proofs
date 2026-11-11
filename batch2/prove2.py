"""Driver: prove/refute Hardin 'Empirical:' recurrences for the families in families.py.
usage: python3 prove2.py OUT.csv [A-number ...]
env: MAXSTATES (raw reachable states cap), MAXWORK (DFS node cap), MAXD (quotient dim cap), PERTURB=1.
"""
import sys, os, re, glob, csv, json, time
from multiprocessing import Pool
sys.setrecursionlimit(10000)
import engine, families

SEQ = '/home/user/corpora/oeisdata/seq'
MAXSTATES = int(os.environ.get('MAXSTATES', '20000'))
MAXWORK = int(os.environ.get('MAXWORK', '4000000'))
MAXD = int(os.environ.get('MAXD', '2500'))
PERTURB = os.environ.get('PERTURB') == '1'

REC = re.compile(r'Empirical:\s*a\(n\)\s*=\s*(.+?)(?:,?\s+for\s+n\s*(>=?)\s*(\d+))?\s*\.?\s*$')
TERM = re.compile(r'([+-]?)\s*(\d*)\s*\*?\s*a\(n-(\d+)\)')
PROOFWORDS = re.compile(r'(?i)\bproo?f|\bproved?\b|confirmed|\bverified|rigorous|\btheorem')

def parse_entry(path):
    d = {'S': '', 'F': [], 'lines': []}
    for line in open(path, encoding='utf-8'):
        tag, rest = line[:2], line[3:].rstrip('\n')
        d['lines'].append(line.rstrip('\n'))
        body = rest.split(' ', 1)[1] if ' ' in rest else ''
        if tag in ('%S', '%T', '%U'): d['S'] += body
        elif tag == '%N': d['N'] = body
        elif tag == '%O': d['O'] = int(body.split(',')[0])
        elif tag == '%F': d['F'].append(body)
    d['terms'] = [int(x) for x in d['S'].split(',') if x.strip()]
    return d

def parse_rec(F):
    for f0 in F:
      for f in (f0, f0.split(';')[0]):
        m = REC.search(f)
        if not m or 'a(n-' not in m.group(1): continue
        body = m.group(1).strip()
        if not body.startswith(('+', '-')): body = '+' + body
        co, used = {}, ''
        for t in TERM.finditer(body):
            j = int(t.group(3))
            co[j] = co.get(j, 0) + int(t.group(2) or 1) * (-1 if t.group(1) == '-' else 1)
            used += t.group(0)
        if re.sub(r'\s', '', used) != re.sub(r'\s', '', body): continue
        n0 = None
        if m.group(3):
            n0 = int(m.group(3)) + (1 if m.group(2) == '>' else 0)  # first n where it is claimed
        return co, n0, f0
    return None

POLY = re.compile(r'Empirical:\s*a\(n\)\s*=\s*([-+*/^()n\d\s]+?)(?:,?\s+for\s+n\s*(>=?)\s*(\d+))?\s*\.?\s*$')
from fractions import Fraction
def parse_poly(F):
    for f0 in F:
        for f in (f0, f0.split(';')[0]):
            m = POLY.search(f)
            if not m or 'n' not in m.group(1): continue
            ex = m.group(1).replace('^', '**')
            if re.search(r'n\s*\(|\)\s*n|\d\s*n|n\s*\d', ex): continue
            py = re.sub(r'\d+', lambda t: 'F(%s)' % t.group(), ex)
            try:
                P = lambda n, py=py: eval(py, {'F': Fraction, '__builtins__': {}}, {'n': Fraction(n)})
                P(3)
            except Exception:
                continue
            exps = [int(e) for e in re.findall(r'n\s*\^\s*(\d+)', m.group(1))]
            deg = max(exps + [1])
            n0 = int(m.group(3)) + (1 if m.group(2) == '>' else 0) if m.group(3) else None
            return P, deg, n0, f0
    return None

def work(A):
    path = f'{SEQ}/{A[:4]}/{A}.seq'
    d = parse_entry(path)
    rec = parse_rec(d['F'])
    base = dict(A=A, family='', states='', order='', formula='', status='', detail='', existing_proof='')
    poly = None
    if rec is None:
        poly = parse_poly(d['F'])
        if poly is None:
            base['status'] = 'no linear recurrence'; return base
        P, deg, n0, ftxt = poly
        co, order = {}, 0      # b(n) = a(n) - P(n) satisfies a recurrence of order D + deg + 1
    else:
        co, n0, ftxt = rec
        order = max(co)
    off = d.get('O', 1)
    start = n0 if n0 is not None else off + order
    hits = [l for l in d['lines'] if l[:2] in ('%F', '%C', '%H', '%Y', '%e') and PROOFWORDS.search(l)]
    base.update(order=order if poly is None else 'polynomial deg %d' % deg, formula=ftxt, existing_proof=('; '.join(h[:120] for h in hits) if hits else 'none'))
    oeis = d['terms']
    last_listed = off + len(oeis) - 1
    tried = []
    for sp in families.specs(d['N']):
        base['family'] = sp['family']
        sh = sp['shift']
        try:
            states, succ, acc = engine.make(sp['w'], sp['m'], sp['ext'], sp['pred'], sp['transpose'], MAXSTATES, MAXWORK)
        except engine.TooLarge as e:
            base.update(status='too large', detail=str(e)); return base
        Tq, fq, c0, D = engine.quotient(succ, acc)
        base['states'] = D
        base['detail'] = f'raw reachable states={len(states)}; quotient D={D}; w={sp["w"]} m={sp["m"]} transpose={sp["transpose"]} shift={sh}'
        if D > MAXD:
            base.update(status='too large', detail=base['detail'] + f'; D>{MAXD}'); return base
        ch_start = max(start, order - sh)          # need n+sh-order >= 0 for the matrix formula
        span = D if poly is None else D + deg + 1  # order of the recurrence satisfied by Q(E)a (resp. a-P)
        nmax = max(last_listed, ch_start + span + 1)
        a = engine.terms(Tq, fq, c0, nmax + sh)    # a[N], N = n + sh
        val = lambda n: a[n + sh]
        mism = [off + i for i, t in enumerate(oeis) if off + i + sh >= 0 and val(off + i) != t]
        if mism:
            tried.append((sp['family'], mism[0], val(mism[0]), oeis[mism[0] - off]))
            continue
        res = {}
        pert = ({**co, max(co): co[max(co)] + 1} if poly is None else 'P+1')
        for label, cof in [('', co)] + ([('PERTURBED ', pert)] if PERTURB else []):
            bad = None
            for n in range(start, ch_start + span + 1):
                if poly is None:
                    rhs = sum(c * val(n - j) for j, c in cof.items())
                else:
                    rhs = P(n) + (1 if cof == 'P+1' else 0)
                if val(n) != rhs:
                    bad = (n, val(n), rhs); break
            res[label] = bad
        bad = res['']
        if bad is None:
            base['status'] = 'PROVED'
            base['detail'] += f'; checked n={start}..{ch_start + span}'
        else:
            base['status'] = f'REFUTED at n={bad[0]}'
            base['detail'] += f'; a({bad[0]})={bad[1]} but recurrence gives {bad[2]}' + (' (n within listed OEIS terms)' if bad[0] <= last_listed else '')
        if PERTURB:
            pb = res['PERTURBED ']
            base['detail'] += '; perturbed: ' + (f'refuted at n={pb[0]}' if pb else 'NOT refuted (!)')
        return base
    if tried:
        base.update(status='TERM MISMATCH', family=tried[0][0],
                    detail='; '.join(f'{f}: first bad n={n} computed {x} listed {y}' for f, n, x, y in tried))
    else:
        base['status'] = 'no spec'
    return base

def candidates():
    done = set(r['A'] for r in csv.DictReader(open('/home/user/oeis-hardin-proofs/proved.csv')))
    out = []
    for p in sorted(glob.glob(SEQ + '/*/*.seq')):
        A = os.path.basename(p)[:-4]
        if A in done: continue
        txt = open(p, encoding='utf-8').read()
        if 'Hardin' not in txt or 'Empirical' not in txt: continue
        N = re.search(r'^%N \S+ (.*)$', txt, re.M)
        if not N or not families.specs(N.group(1)): continue
        out.append(A)
    return out

if __name__ == '__main__':
    OUT = sys.argv[1]
    As = sys.argv[2:] or candidates()
    print(len(As), 'candidates', flush=True)
    t0 = time.time()
    rows = []
    jl = open(OUT + '.jsonl', 'a')
    with Pool(int(os.environ.get('NPROC', '4'))) as pool:
        for r in pool.imap_unordered(work, As, chunksize=1):
            rows.append(r)
            jl.write(json.dumps(r) + '\n'); jl.flush()
            print(r['A'], r['status'], r['family'], r['states'], r['detail'][:150], f'{time.time()-t0:.0f}s', flush=True)
    rows.sort(key=lambda r: r['A'])
    with open(OUT, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=['A', 'family', 'states', 'order', 'formula', 'status', 'existing_proof', 'detail'])
        wr.writeheader(); wr.writerows(rows)
    import collections
    print(collections.Counter(r['status'].split(' at ')[0] for r in rows))
