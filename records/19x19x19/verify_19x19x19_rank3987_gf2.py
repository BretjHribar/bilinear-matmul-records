#!/usr/bin/env python3
"""Sparse exhaustive verifier for decimal-mask GF(2) matmul certificates."""
import argparse, re, sys
from pathlib import Path

def parse(path, dims):
    lines=Path(path).read_text().splitlines(); terms=[]; declared=None
    for line in lines:
        if dims is None:
            m=re.search(r'<(\d+)\s*,\s*(\d+)\s*,\s*(\d+)>',line)
            if m: dims=tuple(map(int,m.groups()))
        f=line.split()
        if not f or f[0].startswith('#'): continue
        if f[0]=='R': f=f[1:]
        if len(f)==1 and f[0].isdigit() and declared is None: declared=int(f[0]);continue
        if len(f)==3 and all(x.isdigit() for x in f):terms.append(tuple(map(int,f)));continue
        raise ValueError(f'bad line: {line[:100]}')
    if dims is None: raise ValueError('dimensions absent; pass --dims N M P')
    if declared is not None and declared != len(terms):raise ValueError(f'declared {declared}, parsed {len(terms)}')
    return dims,terms

def bits(x):
    while x:
        q=x&-x;yield q.bit_length()-1;x^=q

def verify(dims,terms):
    n,m,p=dims; limits=(n*m,m*p,n*p); residual=set(); toggles=0
    seen=set()
    for ti,t in enumerate(terms):
        if t in seen: raise ValueError(f'duplicate term {ti}')
        seen.add(t)
        for x,b in zip(t,limits):
            if not x or x.bit_length()>b:raise ValueError(f'term {ti}: invalid factor')
        ub=list(bits(t[0]));vb=list(bits(t[1]));wb=list(bits(t[2]))
        toggles += len(ub)*len(vb)*len(wb)
        for u in ub:
          for v in vb:
            base=(u*limits[1]+v)*limits[2]
            for w in wb:
              q=base+w
              if q in residual:residual.remove(q)
              else:residual.add(q)
    for i in range(n):
      for j in range(m):
        for k in range(p):
          q=((i*m+j)*limits[1]+(j*p+k))*limits[2]+i*p+k
          if q in residual:residual.remove(q)
          else:residual.add(q)
    return toggles,len(residual)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('certificate');ap.add_argument('--dims',nargs=3,type=int)
    a=ap.parse_args()
    try:
      d,t=parse(a.certificate,tuple(a.dims) if a.dims else None);tog,res=verify(d,t)
    except (OSError,ValueError) as e: print(f'FAIL: {e}',file=sys.stderr);return 1
    print(f'PASS dims=<{d[0]},{d[1]},{d[2]}> rank={len(t)} sparse_toggles={tog} target_coefficients={d[0]*d[1]*d[2]} residual={res}')
    return 0 if res==0 else 1
if __name__=='__main__':raise SystemExit(main())
