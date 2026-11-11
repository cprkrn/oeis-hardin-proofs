"""Parsers: OEIS %N name -> list of candidate specs (one per plausible reading).
A spec = dict(family, w, m, transpose, shift, ext, pred).  N (number of rows of Y) = n + shift.
pred(get) is the anchor predicate; get(di,dj) = X(i+di, j+dj) or None outside the array.
"""
import re, itertools

def dims(name):
    m = re.match(r'Number of (.+?) (?:0\.\.(\d+)|(binary)) (?:arrays|matrices) (.*)$', name)
    if not m:
        return None
    d = m.group(1).replace(' ', '')
    if d.count('X') != 1:
        return None
    a, b = d.split('X')
    def ev(s, n):
        s = s.strip('()')
        if not re.fullmatch(r'[n\d+]+', s):
            raise ValueError
        return eval(s.replace('n', str(n)))
    try:
        if 'n' in a and 'n' not in b:
            transpose, shift, w = False, ev(a, 0), ev(b, 0)
        elif 'n' in b and 'n' not in a:
            transpose, shift, w = True, ev(b, 0), ev(a, 0)
        else:
            return None
    except Exception:
        return None
    mm = 1 if m.group(3) else int(m.group(2))
    return dict(w=w, m=mm, transpose=transpose, shift=shift), m.group(4)

NB = {
    'king-move': [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)],
    'horizontally, diagonally or antidiagonally': [(0, -1), (0, 1), (-1, -1), (1, 1), (-1, 1), (1, -1)],
    'horizontally or antidiagonally': [(0, -1), (0, 1), (-1, 1), (1, -1)],
    'horizontally or vertically': [(0, -1), (0, 1), (-1, 0), (1, 0)],
    'horizontally, vertically or antidiagonally': [(0, -1), (0, 1), (-1, 0), (1, 0), (-1, 1), (1, -1)],
    'horizontally or diagonally': [(0, -1), (0, 1), (-1, -1), (1, 1)],
    'horizontally, vertically or diagonally': [(0, -1), (0, 1), (-1, 0), (1, 0), (-1, -1), (1, 1)],
}
E1 = (-1, 1, -1, 1)

def nums(s):
    return [int(x) for x in re.findall(r'\d+', s)]

def block(get, K, i0=0, j0=0):
    out = []
    for a in range(K):
        row = []
        for b in range(K):
            v = get(i0 + a, j0 + b)
            if v is None:
                return None
            row.append(v)
        out.append(row)
    return out

def matmul(A, B):
    K = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(K)) for j in range(K)] for i in range(K)]

def specs(name):
    r = dims(name)
    if r is None:
        return []
    base, c = r
    out = []
    def add(fam, ext, pred):
        s = dict(base); s.update(family=fam, ext=ext, pred=pred); out.append(s)

    # 1. pattern avoidance, three in a row
    mt = re.fullmatch(r'avoiding (\d \d \d)(?: and (\d \d \d))? horizontally and (\d \d \d)(?: and (\d \d \d))? vertically\.', c)
    if mt:
        H = [tuple(nums(x)) for x in (mt.group(1), mt.group(2)) if x]
        V = [tuple(nums(x)) for x in (mt.group(3), mt.group(4)) if x]
        def pred(g, H=H, V=V):
            h = (g(0, -1), g(0, 0), g(0, 1))
            if None not in h and h in H: return False
            v = (g(-1, 0), g(0, 0), g(1, 0))
            if None not in v and v in V: return False
            return True
        add('avoid 3-patterns h/v', E1, pred)
        return out

    # 2. x equals number of h/v neighbours equal to f(x)
    mt = re.fullmatch(r'with each element x equal to the number its horizontal and vertical neighbors equal to ([\d,]+) for x=([\d,]+)\.', c)
    if mt:
        f = dict(zip(nums(mt.group(2)), nums(mt.group(1))))
        def pred(g, f=f):
            x = g(0, 0)
            t = f[x]
            return x == sum(1 for d in ((0, 1), (0, -1), (1, 0), (-1, 0)) if g(*d) == t)
        add('x = #neighbours equal f(x)', E1, pred)
        return out

    # 3. clockwise perimeter pattern of 3x3 subblocks
    mt = re.fullmatch(r'with each 3 ?X ?3 subblock having clockwise perimeter pattern ([\d or]+)\.', c)
    if mt:
        pats = [tuple(int(ch) for ch in p) for p in re.findall(r'\d+', mt.group(1))]
        if all(len(p) == 8 for p in pats):
            P = [(-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1), (1, 0), (1, -1), (0, -1)]
            rot = set(p[k:] + p[:k] for p in pats for k in range(8))
            for nm, S in (('rotations', rot), ('rotations+reflections', rot | set(q[::-1] for q in rot)), ('fixed start', set(pats))):
                def pred(g, S=S):
                    v = tuple(g(*d) for d in P)
                    return None in v or v in S
                add('3x3 clockwise perimeter (%s)' % nm, E1, pred)
        return out

    # 4. subblock line sums (rows/columns/diagonal/antidiagonal) in / not in sets
    mt = re.match(r'with (?:every|no) (\d) ?[xX] ?(\d) subblock (.*)\.$', c)
    if mt and mt.group(1) == mt.group(2) and ' sum' in mt.group(3):
        K = int(mt.group(1))
        neg0 = c.startswith('with no')
        clauses = re.split(r' and (?=every |no )', mt.group(3))
        rules = []
        okp = True
        for i, cl in enumerate(clauses):
            neg = neg0 if i == 0 else cl.startswith('no ')
            cl = re.sub(r'^(every|no) (\d ?[xX] ?\d )?(subblock )?', '', cl)
            q = re.fullmatch(r'((?:row|column|diagonal|antidiagonal)(?:(?:,| and|, and) (?:row|column|diagonal|antidiagonal))*) sum (not equal to |equal to |not )?([\d ,or]+)', cl)
            if not q:
                okp = False; break
            kinds = re.findall(r'antidiagonal|diagonal|row|column', q.group(1))
            pol = q.group(2) or ''
            if 'not' in pol: neg = True
            rules.append((kinds, set(nums(q.group(3))), neg))
        if okp and rules:
            off = -(K // 2)  # anchor at block centre for odd K, top-left-ish for even
            lo, hi = off, off + K - 1
            def pred(g, K=K, rules=rules, lo=lo):
                B = block(g, K, lo, lo)
                if B is None: return True
                sums = {'row': [sum(r) for r in B], 'column': [sum(B[a][b] for a in range(K)) for b in range(K)],
                        'diagonal': [sum(B[a][a] for a in range(K))], 'antidiagonal': [sum(B[a][K - 1 - a] for a in range(K))]}
                for kinds, S, neg in rules:
                    for k in kinds:
                        for s in sums[k]:
                            if (s in S) == neg: return False
                return True
            add('%dx%d subblock line sums' % (K, K), (lo, hi, lo, hi), pred)
            return out

    # 5. idempotent / commuting KxK subblocks
    mt = re.fullmatch(r'with each (\d) ?X ?(\d) subblock idempotent\.', c)
    if mt:
        K = int(mt.group(1))
        for nm, md in (('integer', None), ('mod 2', 2)):
            def pred(g, K=K, md=md):
                B = block(g, K)
                if B is None: return True
                P = matmul(B, B)
                if md: P = [[x % md for x in r] for r in P]
                return P == B
            add('KxK subblock idempotent (%s)' % nm, (0, K - 1, 0, K - 1), pred)
        return out
    mt = re.fullmatch(r'with every (\d) ?X ?(\d) subblock commuting with each horizontal and vertical neighbor \d ?X ?\d subblock\.', c)
    if mt:
        K = int(mt.group(1))
        def pred(g, K=K):
            A = block(g, K)
            if A is None: return True
            AA = None
            for d in ((0, 1), (1, 0)):
                Bm = block(g, K, *d)
                if Bm is None: continue
                if matmul(A, Bm) != matmul(Bm, A): return False
            return True
        add('KxK subblocks commute with neighbours', (0, K, 0, K), pred)
        return out

    # 6. neighbour-count (equal/unequal), with upper left element zero
    mt = re.fullmatch(r'with every element (equal|unequal) to ([\d, or]+?) (king-move|horizontally[a-z, ]*?) adjacent elements, with upper left element zero\.', c)
    if mt and mt.group(3) in NB:
        eq = mt.group(1) == 'equal'; S = set(nums(mt.group(2))); nb = NB[mt.group(3)]
        def pred(g, eq=eq, S=S, nb=nb):
            x = g(0, 0)
            if g(-1, 0) is None and g(0, -1) is None and x != 0: return False
            cnt = 0
            for d in nb:
                y = g(*d)
                if y is not None and ((y == x) == eq): cnt += 1
            return cnt in S
        add('neighbour count (%s)' % mt.group(3), E1, pred)
        return out

    # 7. each v adjacent to S neighbouring v's
    mt = re.fullmatch(r'with (?:each|every) (\d) (.*?) ?adjacent to ([\d, or]+?) (king-move )?(?:neighboring )?(\d)s\.', c)
    if mt and mt.group(1) == mt.group(5):
        kind = (mt.group(2) or '').strip() or (mt.group(4) or '').strip()
        if kind in NB:
            v = int(mt.group(1)); S = set(nums(mt.group(3))); nb = NB[kind]
            def pred(g, v=v, S=S, nb=nb):
                if g(0, 0) != v: return True
                return sum(1 for d in nb if g(*d) == v) in S
            add('each %s adjacent to S %s' % ('v', kind), E1, pred)
        return out

    # 8. nondecreasing F(pair) in a direction (two clauses)
    FX = r'(x\(i,j\)[-+]x\((?:i-1,j|i,j-1)\)|min\(x\(i,j\),x\((?:i-1,j|i,j-1)\)\)|max\(x\(i,j\),x\((?:i-1,j|i,j-1)\)\)|absolute value of x\(i,j\)-x\((?:i-1,j|i,j-1)\))'
    mt = re.fullmatch(r'with nondecreasing ' + FX + r' in the ([ij]) direction and nondecreasing ' + FX + r' in the ([ij]) direction\.', c)
    if mt:
        def mk(expr, dirc):
            p = (-1, 0) if 'i-1,j' in expr else (0, -1)
            if expr.startswith('min'): f = min
            elif expr.startswith('max'): f = max
            elif expr.startswith('abs'): f = lambda a, b: abs(a - b)
            elif '+' in expr: f = lambda a, b: a + b
            else: f = lambda a, b: a - b
            d = (1, 0) if dirc == 'i' else (0, 1)
            return p, f, d
        R = [mk(mt.group(1), mt.group(2)), mk(mt.group(3), mt.group(4))]
        def pred(g, R=R):
            for p, f, d in R:
                cells = [g(0, 0), g(*p), g(*d), g(d[0] + p[0], d[1] + p[1])]
                if None in cells: continue
                if f(cells[0], cells[1]) > f(cells[2], cells[3]): return False
            return True
        add('nondecreasing pair functions', E1, pred)
        return out
    mt = re.fullmatch(r'with column and row pair sums b\(i,j\)=a\(i,j\)\+a\(i,j-1\) and c\(i,j\)=a\(i,j\)\+a\(i-1,j\) nondecreasing in column and row directions, respectively\.', c)
    if mt:
        def pred(g):
            for p, d in (((0, -1), (1, 0)), ((-1, 0), (0, 1))):
                cells = [g(0, 0), g(*p), g(*d), g(d[0] + p[0], d[1] + p[1])]
                if None in cells: continue
                if cells[0] + cells[1] > cells[2] + cells[3]: return False
            return True
        add('pair sums nondecreasing', E1, pred)
        return out

    # 9. 2x2 subblock predicates
    two = r'(?:2 ?X ?2|2x2)'
    def add2(fam, f):
        def pred(g, f=f):
            a, b, cc, d = g(0, 0), g(0, 1), g(1, 0), g(1, 1)
            if None in (a, b, cc, d): return True
            return f(a, b, cc, d)
        add(fam, (0, 1, 0, 1), pred)
    def diffs(a, b, c_, d):
        return [a - b, a - c_, a - d, b - c_, b - d, c_ - d]
    mt = re.fullmatch(r'with every ' + two + r' subblock having the absolute values of all six edge and diagonal differences no larger than (\d+)\.', c)
    if mt:
        k = int(mt.group(1)); add2('2x2 six |differences| <= k', lambda *v, k=k: max(v) - min(v) <= k); return out
    mt = re.fullmatch(r'with every ' + two + r' subblock having the sum of the absolute values of all six edge and diagonal differences equal to (\d+)\.', c)
    if mt:
        k = int(mt.group(1)); add2('2x2 sum of six |differences| = k', lambda *v, k=k: sum(abs(x) for x in diffs(*v)) == k); return out
    mt = re.fullmatch(r'with every ' + two + r' subblock having the sum of the squares of all six edge and diagonal differences equal to (\d+)\.', c)
    if mt:
        k = int(mt.group(1)); add2('2x2 sum of six squared differences = k', lambda *v, k=k: sum(x * x for x in diffs(*v)) == k); return out
    if re.fullmatch(r'with the minimum plus the upper median equal to the lower median plus the maximum in every ' + two + r' subblock\.', c):
        def f(*v):
            s = sorted(v); return s[0] + s[2] == s[1] + s[3]
        add2('2x2 min+upper median = lower median+max', f); return out
    if re.fullmatch(r'with every ' + two + r' subblock having equal diagonal elements or equal antidiagonal elements\.', c):
        add2('2x2 equal diagonal or antidiagonal', lambda a, b, c_, d: a == d or b == c_); return out
    mt = re.fullmatch(r'with no ' + two + r' subblock having the sum of its diagonal elements greater than the (maximum|absolute difference) of its antidiagonal elements\.', c)
    if mt:
        if mt.group(1) == 'maximum':
            add2('2x2 diag sum <= max antidiag', lambda a, b, c_, d: a + d <= max(b, c_))
        else:
            add2('2x2 diag sum <= |antidiag diff|', lambda a, b, c_, d: a + d <= abs(b - c_))
        return out
    mt = re.fullmatch(r'with every ' + two + r' subblock (diagonal minus antidiagonal sum|sum) nondecreasing (horizontally and vertically|horizontally, vertically and ne-to-sw antidiagonally|horizontally, vertically and antidiagonally ne-to-sw)\.', c)
    if mt:
        if mt.group(1) == 'sum':
            s = lambda g, i, j: None if None in (g(i, j), g(i, j + 1), g(i + 1, j), g(i + 1, j + 1)) else g(i, j) + g(i, j + 1) + g(i + 1, j) + g(i + 1, j + 1)
        else:
            s = lambda g, i, j: None if None in (g(i, j), g(i, j + 1), g(i + 1, j), g(i + 1, j + 1)) else g(i, j) + g(i + 1, j + 1) - g(i, j + 1) - g(i + 1, j)
        pairs = [((0, 0), (0, 1)), ((0, 0), (1, 0))]
        if 'ne-to-sw' in mt.group(2): pairs.append(((0, 1), (1, 0)))
        def pred(g, s=s, pairs=pairs):
            for p, q in pairs:
                x, y = s(g, *p), s(g, *q)
                if x is not None and y is not None and x > y: return False
            return True
        add('2x2 %s nondecreasing %s' % (mt.group(1), mt.group(2)), (0, 2, 0, 2), pred)
        return out
    if re.fullmatch(r'with nondecreasing sum of every two consecutive values in every row and column\.', c):
        def pred(g):
            for d in ((0, 1), (1, 0)):
                a, b = g(0, 0), g(2 * d[0], 2 * d[1])
                if a is not None and b is not None and a > b: return False
            return True
        add('consecutive pair sums nondecreasing', (0, 2, 0, 2), pred)
        return out
    mt = re.fullmatch(r'with the determinants of ' + two + r' subblocks nondecreasing rightwards and downwards\.', c)
    if mt:
        det = lambda g, i, j: None if None in (g(i, j), g(i, j + 1), g(i + 1, j), g(i + 1, j + 1)) else g(i, j) * g(i + 1, j + 1) - g(i, j + 1) * g(i + 1, j)
        def pred(g):
            x = det(g, 0, 0)
            if x is None: return True
            for q in ((0, 1), (1, 0)):
                y = det(g, *q)
                if y is not None and x > y: return False
            return True
        add('2x2 determinants nondecreasing', (0, 2, 0, 2), pred)
        return out
    return out
