import json, csv, os, collections
D = '/tmp/claude-0/-home-user/196a96e7-6caa-565a-8439-ac3cce9abee4/scratchpad/gf/'
S1 = [json.loads(l) for l in open(D + 'stage1.jsonl')]
key = lambda r: (r['A'], r.get('idx', 0))
res = {}
for f in ('s2b.jsonl', 's2c.jsonl', 'stage2a.jsonl'):
    if not os.path.exists(D + f): continue
    for l in open(D + f):
        r = json.loads(l)
        k = (r['A'], r.get('idx', 0))
        if k in res and r['status'] == 'skipped' and res[k]['status'] != 'skipped': continue
        res[k] = r
BRUTE = {
 'A201081': 'listed data confirmed by brute force n=1..9 (3,7,13,25,47,89,169,321,609); g.f. gives 3,5,15,21,...',
 'A194476': 'listed data confirmed by brute force n=1..6; g.f. coefficients are exactly 1/3 of the listed terms (missing factor 3)',
 'A194481': 'listed data confirmed by brute force n=1..6; g.f. coefficients are exactly 1/3 of the listed terms (missing factor 3)',
 'A208410': 'listed data confirmed by brute force n=1..3 (5,187,10682); g.f. gives 2,15,182,...',
 'A232591': 'listed data confirmed by brute force n=1..3 (2,34,12); g.f. coefficients are 1/2 of the listed terms',
 'A280227': 'listed data confirmed by brute force n=1..8; g.f. coefficients are 1/2 of a(n+1) (factor 2 and shift)',
 'A183777': 'a(0)=16 (32 1x5 arrays, halved); g.f. has no constant term. Engine (D=11) PROVES the g.f. for all n>=1 (checked n=1..20); only a(0) is wrong',
}
rows = []
for r in S1:
    k = (r['A'], r['idx'])
    base = dict(A=r['A'], formula=r.get('gf') or r['fline'][:200])
    if r['s1'] == 'UNPARSED':
        rows.append(dict(base, kind='gf', status='skipped', detail='not a parseable rational g.f. (' + r.get('err', '')[:40] + ')')); continue
    if k in res and not (res[k]['status'] == 'skipped' and r['s1'] == 'MISMATCH'):
        x = res[k]; rows.append(dict(base, kind=x['kind'], status=x['status'], detail=x['detail'])); continue
    if r['s1'] == 'MISMATCH':
        sh = r.get('shift_match')
        det = f"g.f. coefficient at n={r['first_bad']} is {r['gfval']} but listed a(n)={r['oeis']}; {r['nbad']}/{r['nterms']} listed terms disagree"
        if sh == -1: det += '; x*G(x) matches all listed terms (g.f. off by factor x / offset convention)'
        elif sh == 1: det += '; G(x)/x matches all listed terms (g.f. has an extra factor x / offset convention)'
        if r['A'] in BRUTE: det += '; ' + BRUTE[r['A']]
        st = f"REFUTED at n={r['first_bad']}" + (' (vs brute-force-confirmed data)' if r['A'] in BRUTE else ' (vs listed data)')
        rows.append(dict(base, kind='gf-vs-data', status=st, detail=det)); continue
    if r['s1'] == 'MATCH_BUT_NONZERO_BELOW_OFFSET':
        rows.append(dict(base, kind='gf-vs-data', status='skipped', detail='matches all listed terms but has nonzero coefficients below the offset (convention issue); no exact model')); continue
    rows.append(dict(base, kind='gf-vs-data', status='skipped', detail=f"no exact model; expansion matches all {r['nterms']} listed terms" + ('; ' + res[k]['detail'][:150] if k in res else '')))
rows.sort(key=lambda r: r['A'])
with open(D + 'gf_results.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['A', 'kind', 'status', 'detail', 'formula']); w.writeheader(); w.writerows(rows)
print(len(rows), collections.Counter((r['kind'].split('(')[0], r['status'].split(' at')[0]) for r in rows))
# Part 2 terms check
trows = []
for l in open(D + 'stage3.jsonl'):
    r = json.loads(l)
    trows.append(dict(A=r['A'], kind='terms-vs-proved-recurrence', status='OK' if not r['bad'] else 'MISMATCH',
                      detail=f"{r['nlisted']} listed terms; n<{r['start']} verified by original prover, rest compared with proved recurrence" + (f"; bad {r['bad'][:5]}" if r['bad'] else ''), formula=''))
seen = set(t['A'] for t in trows)
for x in res.values():
    if x['A'] in seen: continue
    if x['kind'] == 'gf-direct':
        trows.append(dict(A=x['A'], kind='terms-vs-engine', status='OK', detail='all listed terms equal exact transfer-matrix values; ' + x['detail'].split(';')[0], formula='')); seen.add(x['A'])
    elif 'TERM MISMATCH' in x['status']:
        trows.append(dict(A=x['A'], kind='terms-vs-engine', status='MISMATCH', detail=x['detail'], formula='')); seen.add(x['A'])
for A, d in BRUTE.items():
    trows.append(dict(A=A, kind='terms-vs-bruteforce', status='OK', detail=d.split(';')[0], formula=''))
trows.sort(key=lambda r: r['A'])
with open(D + 'terms_check.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=['A', 'kind', 'status', 'detail', 'formula']); w.writeheader(); w.writerows(trows)
print(len(trows), collections.Counter((r['kind'], r['status']) for r in trows))
