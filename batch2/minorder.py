import sys, prove2, engine, families
P = (1 << 61) - 1
def bm(s):
    C, B, L, m, b = [1], [1], 0, 1, 1
    for n in range(len(s)):
        d = s[n] % P
        for i in range(1, L + 1): d = (d + C[i] * s[n - i]) % P
        if d == 0: m += 1; continue
        T = C[:]; coef = d * pow(b, P - 2, P) % P
        C = C + [0] * (len(B) + m - len(C))
        for i in range(len(B)): C[i + m] = (C[i + m] - coef * B[i]) % P
        if 2 * L <= n: L, B, b, m = n + 1 - L, T, d, 1
        else: m += 1
    return L
for A in sys.argv[1:]:
    d = prove2.parse_entry(f'{prove2.SEQ}/{A[:4]}/{A}.seq')
    sp = [s for s in families.specs(d['N'])][0]
    st, succ, acc = engine.make(sp['w'], sp['m'], sp['ext'], sp['pred'], sp['transpose'], 10**6, 10**9)
    Tq, fq, c0, D = engine.quotient(succ, acc)
    a = engine.terms(Tq, fq, c0, 2 * D + 40)
    print(A, 'D', D, 'BM order mod p (from N=1)', bm(a[1:]), 'from N=60', bm(a[60:]))
