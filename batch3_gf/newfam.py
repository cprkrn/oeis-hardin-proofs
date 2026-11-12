"""Additional local-constraint families (same spec format as families.py).
A spec may carry 'preds': a list of predicates; the count is then the SUM of the counts for each predicate
(used for disjoint unions, e.g. 'all 2x2 subblock sums the same' = sum over the common value s)."""
import re
from families import dims, block, matmul

def det(B):
    K = len(B)
    if K == 1: return B[0][0]
    if K == 2: return B[0][0] * B[1][1] - B[0][1] * B[1][0]
    return sum((-1) ** j * B[0][j] * det([r[:j] + r[j + 1:] for r in B[1:]]) for j in range(K))

def specs(name):
    out = []
    # 1D: length n+k 0..m arrays with at most one/two downsteps in every K consecutive neighbor pairs
    m1 = re.fullmatch(r'Number of length n\+(\d+) 0\.\.(\d+) arrays with at most (one|two) downsteps? in every (\d+) consecutive neighbor pairs\.', name)
    if m1:
        k, m, lim, K = int(m1.group(1)), int(m1.group(2)), {'one': 1, 'two': 2}[m1.group(3)], int(m1.group(4))
        def pred(g, K=K, lim=lim):
            v = [g(i, 0) for i in range(K + 1)]
            if None in v: return True
            return sum(v[i + 1] < v[i] for i in range(K)) <= lim
        out.append(dict(family='1D downsteps per window', w=1, m=m, transpose=False, shift=k, ext=(0, K, 0, 0), pred=pred))
        return out
    r = dims(name)
    if r is None: return out
    base, c = r
    def add(fam, ext, pred=None, preds=None):
        s = dict(base); s.update(family=fam, ext=ext, pred=pred)
        if preds: s['preds'] = preds
        out.append(s)
    mt = re.fullmatch(r'with (?:every|each) (\d) X (\d) subblock (singular|nonsingular)\.', c)
    if mt:
        K = int(mt.group(1)); sing = mt.group(3) == 'singular'
        def pred(g, K=K, sing=sing):
            B = block(g, K)
            if B is None: return True
            return (det(B) == 0) == sing
        add('KxK subblock (non)singular', (0, K - 1, 0, K - 1), pred); return out
    mt = re.fullmatch(r'with (every|no) (\d) X (\d) subblock (having nonzero determinant and )?commuting with (?:each of its|any of its|every) horizontal (?:and|or) vertical (?:\d X \d subblock )?neighbors?\.', c)
    if mt:
        every = mt.group(1) == 'every'; K = int(mt.group(2)); nz = bool(mt.group(4))
        for stepname, st in (('overlapping', 1), ('adjacent', K)):
            def pred(g, K=K, st=st, every=every, nz=nz):
                A = block(g, K)
                if A is None: return True
                if nz and det(A) == 0: return False
                for d in ((0, st), (st, 0), (0, -st), (-st, 0)):
                    Bm = block(g, K, *d)
                    if Bm is None: continue
                    com = matmul(A, Bm) == matmul(Bm, A)
                    if every and not com: return False
                    if not every and com: return False
                return True
            lo, hi = -st, K - 1 + st
            add('KxK commute (%s neighbours)' % stepname, (lo, hi, lo, hi), pred)
        return out
    mt = re.fullmatch(r'avoiding patterns (\d{3}) and (\d{3}) in rows(?:,| and) columns( and nw-to-se diagonals)?\.', c)
    if mt:
        P = {tuple(int(ch) for ch in mt.group(1)), tuple(int(ch) for ch in mt.group(2))}
        diag = bool(mt.group(3))
        def pred(g, P=P, diag=diag):
            dirs = [(0, 1), (1, 0)] + ([(1, 1)] if diag else [])
            for a, b in dirs:
                t = (g(-a, -b), g(0, 0), g(a, b))
                if None not in t and t in P: return False
            return True
        add('avoid two 3-patterns rows/cols' + ('/diag' if diag else ''), (-1, 1, -1, 1), pred); return out
    mt = re.fullmatch(r'with all 2 X 2 subblock sums the same\.', c)
    if mt:
        preds = []
        for s in range(4 * base['m'] + 1):
            def pred(g, s=s):
                B = block(g, 2)
                return B is None or sum(map(sum, B)) == s
            preds.append(pred)
        add('all 2x2 sums equal (sum over s)', (0, 1, 0, 1), preds=preds); return out
    return out
