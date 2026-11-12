"""Parse an OEIS 'Empirical g.f.:' text into an exact rational function N(x)/Q(x) (integer/Fraction coeffs)."""
import re
from fractions import Fraction

def padd(a, b):
    n = max(len(a), len(b)); return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])
def pneg(a): return [-c for c in a]
def pmul(a, b):
    if not a or not b: return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): r[i + j] += x * y
    return trim(r)
def trim(a):
    a = list(a)
    while a and a[-1] == 0: a.pop()
    return a

class R:  # rational function
    def __init__(s, n, d=None): s.n, s.d = trim(n), trim(d if d is not None else [1])
    def __add__(s, o): return R(padd(pmul(s.n, o.d), pmul(o.n, s.d)), pmul(s.d, o.d))
    def __sub__(s, o): return R(padd(pmul(s.n, o.d), pneg(pmul(o.n, s.d))), pmul(s.d, o.d))
    def __mul__(s, o): return R(pmul(s.n, o.n), pmul(s.d, o.d))
    def __truediv__(s, o):
        if not o.n: raise ZeroDivisionError
        return R(pmul(s.n, o.d), pmul(s.d, o.n))
    def __pow__(s, k):
        r = R([1])
        for _ in range(k): r = r * s
        return r

TOK = re.compile(r'\s*(\d+|x|[-+*/^()])')

def parse(text):
    s = text.strip()
    pos, toks = 0, []
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m:
            if s[pos:].strip() == '': break
            raise ValueError('bad char at %r' % s[pos:pos + 10])
        toks.append(m.group(1)); pos = m.end()
    i = [0]
    def peek(): return toks[i[0]] if i[0] < len(toks) else None
    def eat(t=None):
        tk = peek()
        if t is not None and tk != t: raise ValueError('expected %s got %s' % (t, tk))
        i[0] += 1; return tk
    def expr():
        if peek() in ('+', '-'):
            sg = eat(); v = term(); v = R(pneg(v.n), v.d) if sg == '-' else v
        else: v = term()
        while peek() in ('+', '-'):
            op = eat(); t = term(); v = v + t if op == '+' else v - t
        return v
    def term():
        v = factor()
        while True:
            if peek() in ('*', '/'):
                op = eat(); f = factor(); v = v * f if op == '*' else v / f
            elif peek() in ('(', 'x') or (peek() and peek().isdigit()):
                v = v * factor()   # implicit multiplication
            else: return v
    def factor():
        b = atom()
        if peek() == '^':
            eat()
            if peek() == '(':
                eat(); e = int(eat()); eat(')')
            else: e = int(eat())
            b = b ** e
        return b
    def atom():
        tk = eat()
        if tk == '(':
            v = expr(); eat(')'); return v
        if tk == 'x': return R([0, 1])
        if tk and tk.isdigit(): return R([int(tk)])
        if tk in ('-',):
            v = factor(); return R(pneg(v.n), v.d)
        raise ValueError('unexpected %r' % tk)
    v = expr()
    if i[0] != len(toks): raise ValueError('trailing tokens')
    return v

def series(r, nterms):
    """Power series coefficients c_0..c_{nterms-1} of N/Q; requires Q(0) != 0."""
    N, Q = r.n, r.d
    if not Q or Q[0] == 0: raise ValueError('Q(0)=0')
    q0 = Fraction(Q[0]); c = []
    for n in range(nterms):
        v = Fraction(N[n] if n < len(N) else 0)
        for j in range(1, min(n, len(Q) - 1) + 1): v -= Q[j] * c[n - j]
        c.append(v / q0)
    return c

GFRE = re.compile(r'Empirical g\.f\.(?: for [^:]*)?:?\s*(?:G\.f\.:\s*)?(.*)$')

def extract(fline):
    """Return the g.f. expression text from a %F body containing 'Empirical g.f.:' or None."""
    m = GFRE.search(fline)
    if not m: return None
    t = m.group(1)
    t = re.split(r'\s*\.?\s*-\s*_|\.\s+_|\s*\[|\s+\(End\)|;|,\s*which| with | and | = |\s\.\s|\.\s+a\(n\)', t)[0]
    t = t.strip().rstrip('.').strip()
    return t
