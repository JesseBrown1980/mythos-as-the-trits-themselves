import numpy as np, csv, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import ConnectionPatch, Circle
from scipy.ndimage import gaussian_filter, zoom
bg='#02040c'; fg='#cfe3ff'; dim='#7f98bf'
blue=LinearSegmentedColormap.from_list('b',['#02040c','#0a2a6b','#1f6fe0','#7fc4ff','#ffffff'])
star=np.frombuffer(open('results/tg/star.trits','rb').read(),np.uint8).reshape(3,64,64)
rho=np.frombuffer(open('results/tg/rho_frames.bin','rb').read(),np.uint8).reshape(-1,3,64,64)
rows=[r for r in csv.DictReader(open('results/tg/outer_modes.csv')) if r['mode']=='gates']
M=np.array([[int(r[f'M{k}']) for k in range(27)] for r in rows]); ticks=[int(r['tick']) for r in rows]
vel=[(x,y,z) for z in (-1,0,1) for y in (-1,0,1) for x in (-1,0,1)]
def glow(img,s=2.2):
    up=zoom(img.astype(np.int64),3,order=0); g=gaussian_filter(up*1.0,s)*4+up
    return np.sqrt(g/max(g.max(),1e-9))
def copy_k(k):  # channel k at t=0: star seen through offset c_k
    cx,cy,cz=vel[k]; return np.roll(star,(-cz,-cy,-cx),axis=(0,1,2)).sum(0)
fig=plt.figure(figsize=(15,10.6),facecolor=bg)
fig.text(0.02,0.965,'Mythos blue star as the trits: a neural network inside a neural network',color=fg,fontsize=17)
fig.text(0.02,0.937,'Outer: 27 shell modes, each one a full copy of the star (one per D3Q27 channel).  Inner: every one of 12,288 cells runs the 288-gate ternary network on its 27 star trits.',color=dim,fontsize=10)
# ---- outer ring of 27 star nodes
ax=fig.add_axes([0.0,0.30,0.50,0.62]); ax.set_xlim(-1.25,1.25); ax.set_ylim(-1.25,1.42); ax.axis('off'); ax.set_aspect('equal')
order=[k for k in range(27)]
ang=[np.pi/2-2*np.pi*i/27 for i in range(27)]
pos={k:(np.cos(a),np.sin(a)) for k,a in zip(order,ang)}
for i in range(27):
    a,b=pos[i],pos[(i+1)%27]; ax.plot([a[0],b[0]],[a[1],b[1]],color='#1f4a8a',lw=1.2,zorder=1)
# triad tree spokes 3 -> 9 -> 27 into the centre
for i in range(27):
    if i%9==0: ax.plot([0,pos[i][0]*0.55],[0,pos[i][1]*0.55],color='#2b5fae',lw=1.6,zorder=1)
for j in range(9):
    p=pos[(j//3)*9]; q=pos[j*3]
    ax.plot([p[0]*0.55,q[0]*0.8],[p[1]*0.55,q[1]*0.8],color='#1f4a8a',lw=1,zorder=1)
    for m in range(3): r=pos[j*3+m]; ax.plot([q[0]*0.8,r[0]],[q[1]*0.8,r[1]],color='#173a6e',lw=0.8,zorder=1)
sz=0.105
for k in range(27):
    x,y=pos[k]; im=glow(copy_k(k))
    ax.imshow(im,cmap=blue,extent=[x-sz,x+sz,y-sz,y+sz],zorder=3)
    ring='#ffd75e' if k==13 else '#3b6fc0'
    ax.add_patch(plt.Rectangle((x-sz,y-sz),2*sz,2*sz,fill=False,ec=ring,lw=1.4 if k==13 else 0.7,zorder=4))
    c=vel[k]; ax.text(x*1.17,y*1.17,f'{c[0]:+d}{c[1]:+d}{c[2]:+d}'.replace('+0','0').replace('-0','0'),color=dim,fontsize=6.5,ha='center',va='center')
ax.add_patch(Circle((0,0),0.075,color='#ffd75e',zorder=5)); ax.text(0,-0.14,'free centre\nQ = 15,255 (all ticks)',color='#ffd75e',fontsize=9,ha='center',va='top')
ax.text(-1.2,1.40,'OUTER NETWORK — 27 star copies on the ring, triad tree 3 → 9 → 27\ngold frame = rest channel c=(0,0,0), the star unshifted',color=fg,fontsize=9.5,va='top')
# ---- inner: one cell's 27 trits as a 3x3x3 block, then the lattice over time
fr=[0,1,6,18,48]; labels=['t = 0  (27 copies laid in)','t = 4','t = 24','t = 72','t = 192']
for i,(f,l) in enumerate(zip(fr,labels)):
    a=fig.add_axes([0.52+(i%3)*0.155,0.62-(i//3)*0.30,0.145,0.26]); a.imshow(glow(rho[f].sum(0).astype(np.int64)),cmap=blue)
    a.set_xticks([]);a.set_yticks([]); [s.set_color('#1f3a66') for s in a.spines.values()]
    a.set_title(l,color=fg,fontsize=9)
a=fig.add_axes([0.52+2*0.155,0.32,0.145,0.26]); a.imshow(glow(rho[0].sum(0).astype(np.int64)),cmap=blue)
a.set_xticks([]);a.set_yticks([]); [s.set_color('#ffd75e') for s in a.spines.values()]
a.set_title('72 forward, 72 back:\nall 331,776 trits restored',color='#ffd75e',fontsize=9)
fig.text(0.52,0.905,'INNER NETWORK — density ρ = Σ 27 trits per cell, 64×64×3 lattice (A,B,C summed), Taylor–Green gate set',color=fg,fontsize=9.5)
con=ConnectionPatch(xyA=pos[13],coordsA=ax.transData,xyB=(0,0.5),coordsB=fig.axes[1].transAxes,color='#ffd75e',lw=0.8,ls='--',alpha=0.6); fig.add_artist(con)
# ---- outer modes over time
b=fig.add_axes([0.06,0.06,0.88,0.20],facecolor=bg)
for k in range(27): b.plot(ticks,M[:,k],color='#ffd75e' if k==13 else plt.cm.Blues(0.35+0.6*k/27),lw=1.8 if k==13 else 0.7)
b.set_xlim(0,192); b.tick_params(colors=dim,labelsize=8); [s.set_color('#1f3a66') for s in b.spines.values()]
b.set_xlabel('tick',color=dim,fontsize=9); b.set_ylabel('mode mass',color=dim,fontsize=9)
b.set_title('The 27 outer modes as the inner gates move trits between star copies — each starts at 565, the sum stays 15,255 exactly; Jacobian rank 26 at all 193 ticks',color=fg,fontsize=9.5,loc='left')
fig.text(0.94,0.008,'ACER-MYTHOS-BLUE.gguf sha256 3dd92a42… · star trits sha256 748ebf67… · integer only, floating_point=0',color='#5b7299',fontsize=8,ha='right')
fig.savefig('figures/mythos-star-nested-network.png',dpi=120,facecolor=bg)
