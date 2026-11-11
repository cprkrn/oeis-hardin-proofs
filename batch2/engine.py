"""Generic transfer-matrix engine for fixed-width arrays with a local constraint.

The true array X has one fixed dimension (width w) and one varying dimension (N).
We generate Y = X (rows vary) or Y = X^T (columns vary) row by row.
A constraint is a family of 'anchor' predicates, one per cell (i,j) of X, each a function of
the cells X(i+di, j+dj) with di in [r0,r1], dj in [c0,c1]; out-of-range cells read as None
(the predicate itself decides how to treat the boundary).  The array is valid iff every
anchor predicate is true.

State (explicit definition): the last L = up+down rows of Y, padded with None above row 1,
where up/down are the Y-row extents of the anchor neighbourhood.  When a new row is appended,
all anchors in the row 'down' rows above it are checked (their neighbourhood is now complete).
At the end 'down' None rows are appended and checked.  Hence for every N >= 0
    a(N) = e_{s0}^T T^N f,
T = transition matrix on reachable states, f = acceptance vector.
We then compute the coarsest partition of the reachable states that is stable for
(acceptance, number of transitions into each block) (a counting bisimulation).  The span of
block-indicator vectors contains f and is T-invariant, so a(N) = e^T Tq^N fq exactly, with Tq
the D x D quotient matrix (D = number of blocks <= number of reachable states).  By Cayley-
Hamilton, a(N) obeys a linear recurrence of order <= D valid for all N >= 0.
"""
import collections


class TooLarge(Exception):
    pass


def make(w, m, ext, pred, transpose, maxstates=20000, maxwork=3_000_000):
    r0, r1, c0, c1 = ext
    if transpose:
        yr0, yr1, yc0, yc1 = c0, c1, r0, r1
    else:
        yr0, yr1, yc0, yc1 = r0, r1, c0, c1
    up, down = max(0, -yr0), max(0, yr1)
    left, right = max(0, -yc0), max(0, yc1)
    L = up + down
    vals = range(m + 1)

    def check_anchor(window, c):
        arow = window[up]
        if transpose:
            def get(di, dj):
                rr = window[up + dj] if 0 <= up + dj < len(window) else None
                cc = c + di
                if rr is None or cc < 0 or cc >= w:
                    return None
                return rr[cc]
        else:
            def get(di, dj):
                rr = window[up + di] if 0 <= up + di < len(window) else None
                cc = c + dj
                if rr is None or cc < 0 or cc >= w:
                    return None
                return rr[cc]
        return pred(get)

    work = [0]

    def successors(state):
        out = []
        new = []
        window = list(state) + [new]
        anchor_real = window[up] is not None

        def rec(cc):
            if cc == w:
                out.append(tuple(new))
                return
            for v in vals:
                work[0] += 1
                new.append(v)
                ok = True
                if anchor_real:
                    cs = [cc - right] if cc < w - 1 else range(cc - right, w)
                    for c in cs:
                        if c >= 0 and not check_anchor(window, c):
                            ok = False
                            break
                if ok:
                    rec(cc + 1)
                new.pop()
        rec(0)
        return out

    def accept(state):
        st = list(state)
        for _ in range(down):
            window = st + [None]
            if window[up] is not None:
                for c in range(w):
                    if not check_anchor(window, c):
                        return False
            st = window[1:]
        return True

    s0 = (None,) * L
    idx = {s0: 0}
    states = [s0]
    succ = []
    i = 0
    while i < len(states):
        st = states[i]
        lst = []
        for row in successors(st):
            ns = (st + (row,))[1:] if L else ()
            j = idx.get(ns)
            if j is None:
                j = idx[ns] = len(states)
                states.append(ns)
                if len(states) > maxstates:
                    raise TooLarge('states>%d' % maxstates)
            lst.append(j)
        succ.append(lst)
        if work[0] > maxwork:
            raise TooLarge('work>%d' % maxwork)
        i += 1
    acc = [1 if accept(s) else 0 for s in states]
    return states, succ, acc


def bisim(succ, acc):
    n = len(succ)
    cls = list(acc)
    ncls = len(set(cls))
    while True:
        sig = {}
        new = [0] * n
        for s in range(n):
            key = (cls[s], tuple(sorted(collections.Counter(cls[t] for t in succ[s]).items())))
            new[s] = sig.setdefault(key, len(sig))
        if len(sig) == ncls:
            return new, ncls
        cls, ncls = new, len(sig)


def quotient(succ, acc, maxD=None):
    cls, D = bisim(succ, acc)
    rep = {}
    for s, c in enumerate(cls):
        rep.setdefault(c, s)
    Tq = [None] * D
    fq = [0] * D
    for c, s in rep.items():
        Tq[c] = list(collections.Counter(cls[t] for t in succ[s]).items())
        fq[c] = acc[s]
    return Tq, fq, cls[0], D


def terms(Tq, fq, c0, Nmax):
    """a(N) for N = 0..Nmax via g_{N+1} = Tq g_N, g_0 = fq; a(N) = g_N[c0]."""
    g = list(fq)
    out = [g[c0]]
    for _ in range(Nmax):
        g = [sum(k * g[t] for t, k in row) for row in Tq]
        out.append(g[c0])
    return out


def brute(Nrows, w, m, ext, pred, transpose):
    """Direct count of valid arrays (independent check for tiny sizes)."""
    import itertools
    R, C = (w, Nrows) if transpose else (Nrows, w)
    cnt = 0
    for flat in itertools.product(range(m + 1), repeat=R * C):
        X = [flat[i * C:(i + 1) * C] for i in range(R)]
        ok = True
        for i in range(R):
            for j in range(C):
                def get(di, dj, i=i, j=j):
                    a, b = i + di, j + dj
                    return X[a][b] if 0 <= a < R and 0 <= b < C else None
                if not pred(get):
                    ok = False
                    break
            if not ok:
                break
        cnt += ok
    return cnt
