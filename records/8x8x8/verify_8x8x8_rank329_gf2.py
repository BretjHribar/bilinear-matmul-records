#!/usr/bin/env python3
"""Independent UV-slice verifier for the GF(2) <8,8,8> rank-329 scheme."""
import hashlib,json,sys
from pathlib import Path
N=M=P=8; EXPECTED=329; AB=N*M; BC=M*P; AC=N*P
def bits(x):
    while x:
        q=x&-x;yield q.bit_length()-1;x^=q
def parse(path):
    rank=None;terms=[]
    for raw in Path(path).read_text().splitlines():
        f=raw.split()
        if not f or f[0].startswith('#'):continue
        if f[0]=='R':f=f[1:]
        if rank is None and len(f)==1:rank=int(f[0]);continue
        if len(f)!=3:raise ValueError(f'bad certificate line: {raw[:80]}')
        terms.append(tuple(map(int,f)))
    return rank,terms
def main():
    path=Path(sys.argv[1] if len(sys.argv)>1 else 'matmul_8x8x8_rank329_gf2.txt')
    rank,terms=parse(path); residual=[0]*(AB*BC); seen=set();expanded=0
    valid=rank==EXPECTED and len(terms)==EXPECTED
    for t in terms:
        valid &= t not in seen and all(x>0 and x.bit_length()<=b for x,b in zip(t,(AB,BC,AC)));seen.add(t)
        us=list(bits(t[0]));vs=list(bits(t[1]));expanded+=len(us)*len(vs)
        for u in us:
            for v in vs:residual[u*BC+v]^=t[2]
    for i in range(N):
      for j in range(M):
        for k in range(P):residual[(i*M+j)*BC+j*P+k]^=1<<(i*P+k)
    rc=sum(x.bit_count() for x in residual);valid &= rc==0
    print(json.dumps({'tensor':[N,M,P],'rank':len(terms),'declared_rank':rank,'expanded_uv_pairs':expanded,
      'target_coefficients':N*M*P,'brent_coefficients':AB*BC*AC,'residual_coefficients':rc,
      'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'verified':bool(valid)},indent=2,sort_keys=True))
    return 0 if valid else 1
if __name__=='__main__':raise SystemExit(main())
