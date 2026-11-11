import sys,re,itertools
sys.path.insert(0,'/home/user/oeis-hardin-proofs')
def parse(A):
    t=open(f'/home/user/corpora/oeisdata/seq/{A[:4]}/{A}.seq',encoding='utf-8').read()
    S=''.join(l.split(' ',2)[2] for l in t.splitlines() if l[:2] in('%S','%T','%U'))
    terms=[int(x) for x in S.split(',') if x]
    F=[l.split(' ',2)[2] for l in t.splitlines() if l.startswith('%F') and 'Empirical' in l][0]
    body=F.split('=',1)[1]
    co={int(j):int((s+(c or '1')).replace(' ','')) for s,c,j in re.findall(r'([+-]?)\s*(\d*)\*?a\(n-(\d+)\)',' +'+body.strip())}
    return terms,co
def count(k,f,nmax):
    vals=range(5); rows=list(itertools.product(vals,repeat=k))
    def ok(a,b,c):
        for j in range(k):
            x=b[j]; tgt=f[x]; cnt=0
            for nb in ([a[j]] if a else [])+([c[j]] if c else [])+([b[j-1]] if j>0 else [])+([b[j+1]] if j<k-1 else []):
                cnt+= nb==tgt
            if cnt!=x: return False
        return True
    res={1:sum(ok(None,r,None) for r in rows)}
    vec={}
    for a in rows:
        for b in rows:
            if ok(None,a,b): vec[(a,b)]=vec.get((a,b),0)+1
    cache={}
    for n in range(2,nmax+1):
        res[n]=sum(c for (a,b),c in vec.items() if ok(a,b,None))
        nv={}
        for (a,b),c in vec.items():
            if (a,b) not in cache: cache[(a,b)]=[x for x in rows if ok(a,b,x)]
            for x in cache[(a,b)]: nv[(b,x)]=nv.get((b,x),0)+c
        vec=nv
    return res
for A,k,f,n0 in [('A197230',3,[4,3,0,1,2],23),('A196074',4,[0,3,2,1,4],38)]:
    terms,co=parse(A); order=max(co)
    a=count(k,f,n0+3)
    ok=all(a[i+1]==t for i,t in enumerate(terms))
    print(A,'order',order,'OEIS terms reproduced:',ok,len(terms))
    for n in range(n0,n0+3):
        pred=sum(c*a[n-j] for j,c in co.items())
        print(f'  n={n} actual={a[n]} recurrence={pred} {"OK" if pred==a[n] else "FAILS"}')
