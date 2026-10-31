import struct,sys,json
import numpy as np
p=sys.argv[1]; b=open(p,'rb').read(); o=0
def rd(f):
    global o; v=struct.unpack_from(f,b,o); o+=struct.calcsize(f); return v[0]
def s():
    global o; n=rd('<Q'); v=b[o:o+n].decode('utf8','replace'); o+=n; return v
def val(t):
    if t in (0,1): return rd('<B' if t==0 else '<b')
    if t in (2,3): return rd('<H' if t==2 else '<h')
    if t in (4,5): return rd('<I' if t==4 else '<i')
    if t==6: return rd('<f')
    if t==7: return bool(rd('<B'))
    if t==8: return s()
    if t==9:
        et=rd('<I'); n=rd('<Q'); return [val(et) for _ in range(n)]
    if t in (10,11): return rd('<Q' if t==10 else '<q')
    if t==12: return rd('<d')
assert b[:4]==b'GGUF'; o=4
ver=rd('<I'); nt=rd('<Q'); nk=rd('<Q')
kv={}
for _ in range(nk):
    k=s(); t=rd('<I'); kv[k]=val(t)
ti=[]
for _ in range(nt):
    n=s(); nd=rd('<I'); dims=[rd('<Q') for _ in range(nd)]; t=rd('<I'); off=rd('<Q'); ti.append((n,dims,t,off))
al=kv.get('general.alignment',32); base=(o+al-1)//al*al
print('version',ver,'tensors',nt,'kv',nk,'data_base',base)
for k,v in kv.items(): print(' ',k,'=',v)
T={}
for n,dims,t,off in ti:
    cnt=int(np.prod(dims)); a=np.frombuffer(b,dtype='<f4',count=cnt,offset=base+off)
    T[n]=a.reshape(dims[::-1]) if len(dims)>1 else a
    print(f'  {n} dims={dims} type={t} off={off} min={a.min():.6g} max={a.max():.6g} mean={a.mean():.6g} nonzero={int((a!=0).sum())}/{cnt} uniq={len(np.unique(a))}')
np.savez(sys.argv[2] if len(sys.argv)>2 else 'star.npz',**T)
print('end of data', base+max(off+4*int(np.prod(d)) for _,d,_,off in ti), 'file', len(b))
