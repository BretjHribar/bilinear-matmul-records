#!/usr/bin/env python3
from pathlib import Path

def load(path):
    out=[]
    for line in Path(path).read_text().splitlines():
        f=line.split()
        if f and f[0]=='R':f=f[1:]
        if len(f)==3 and all(x.isdigit() for x in f):out.append(tuple(map(int,f)))
    return out

def kron(a,ar,ac,b,br,bc):
    y=0
    for i in range(ar):
      for j in range(ac):
       if a>>(i*ac+j)&1:
        for k in range(br):
         for l in range(bc):
          if b>>(k*bc+l)&1:y ^= 1<<((i*br+k)*(ac*bc)+j*bc+l)
    return y

root=Path(__file__).resolve().parent/'sources'
a=load(root/'matmul_4x4_rank47_d450_gf2.txt')
b=load(root/'matmul_2x2_rank7_strassen_gf2.txt')
terms=[(kron(u1,4,4,u2,2,2),kron(v1,4,4,v2,2,2),kron(w1,4,4,w2,2,2))
       for u1,v1,w1 in a for u2,v2,w2 in b]
assert len(a)==47 and len(b)==7 and len(terms)==329
with open(Path(__file__).resolve().parent/'matmul_8x8x8_rank329_gf2.txt','w') as f:
    f.write('# Exact bilinear matrix-multiplication certificate over GF(2)\n')
    f.write('# tensor <8,8,8>: (8x8) times (8x8) -> (8x8)\n')
    f.write('# rank 329 = 47 * 7 field-aware tensor composition\n329\n')
    for u,v,w in terms:f.write(f'{u} {v} {w}\n')
