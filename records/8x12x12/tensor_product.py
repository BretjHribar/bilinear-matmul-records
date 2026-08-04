#!/usr/bin/env python3
"""Materialize a Kronecker/tensor product of two decimal-mask schemes."""
from pathlib import Path

def load(path):
    terms=[]
    for line in Path(path).read_text().splitlines():
        f=line.split()
        if f and f[0]=='R':f=f[1:]
        if len(f)==3 and all(x.isdigit() for x in f):terms.append(tuple(map(int,f)))
    return terms

def kron_mask(a,ar,ac,b,br,bc):
    out=0
    for i in range(ar):
      for j in range(ac):
        if a>>(i*ac+j)&1:
          for k in range(br):
            for l in range(bc):
              if b>>(k*bc+l)&1:
                out ^= 1 << (((i*br+k)*(ac*bc)) + (j*bc+l))
    return out

def product(left,right,d1,d2):
    n1,m1,p1=d1;n2,m2,p2=d2
    for u1,v1,w1 in left:
      for u2,v2,w2 in right:
        yield (kron_mask(u1,n1,m1,u2,n2,m2),
               kron_mask(v1,m1,p1,v2,m2,p2),
               kron_mask(w1,n1,p1,w2,n2,p2))

if __name__=='__main__':
    root=Path(__file__).resolve().parent/'sources'
    left=load(root/'matmul_4x4_rank47_d450_gf2.txt')
    right=load(root/'matmul_2x3x3_rank15_catalog_gf2.txt')
    terms=list(product(left,right,(4,4,4),(2,3,3)))
    assert len(left)==47 and len(right)==15 and len(terms)==705
    out=Path(__file__).resolve().parent/'matmul_8x12x12_rank705_gf2.txt'
    with out.open('w') as f:
      f.write('# Exact bilinear matrix-multiplication certificate over GF(2)\n')
      f.write('# tensor <8,12,12>: (8x12) times (12x12) -> (8x12)\n')
      f.write('# rank 705 = 47 * 15 tensor composition\n705\n')
      for u,v,w in terms:f.write(f'{u} {v} {w}\n')
