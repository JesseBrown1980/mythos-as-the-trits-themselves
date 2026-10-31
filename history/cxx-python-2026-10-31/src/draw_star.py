import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage import gaussian_filter, zoom
import sys; T=np.load(sys.argv[1]); S=[T['slice_'+k] for k in 'ABC']; R=[T['radial_'+k] for k in 'ABC']
blue=LinearSegmentedColormap.from_list('b',['#02040c','#0a2a6b','#1f6fe0','#7fc4ff','#ffffff'])
bg='#02040c'
fig=plt.figure(figsize=(13,9.6),facecolor=bg)
# --- the star: all three slices summed, upsampled, with glow
tot=sum(S); up=zoom(tot,8,order=0)
glow=gaussian_filter(up,6)*6+gaussian_filter(up,22)*12+gaussian_filter(up,1.2)
img=np.sqrt(glow/glow.max())
ax=fig.add_axes([0.03,0.36,0.55,0.58]); ax.imshow(img,cmap=blue,origin='upper'); ax.axis('off')
cy,cx=32*8+4,32*8+4; ax.plot(cx,cy,'+',color='#ffd75e',ms=14,mew=1.2)
ax.set_title('ACER-MYTHOS-BLUE  ·  slices A+B+C, centre [32,32]',color='#cfe3ff',fontsize=12,loc='left')
# --- three-colour view: A, B, C on separate channels
ax2=fig.add_axes([0.60,0.36,0.38,0.58])
rgb=np.stack([np.sqrt(gaussian_filter(zoom(s,8,order=0),3)) for s in S],-1)
rgb/=rgb.max(); ax2.imshow(rgb,origin='upper'); ax2.axis('off')
ax2.set_title('A = red, B = green, C = blue',color='#cfe3ff',fontsize=11,loc='left')
# --- the three raw 64x64 slices
for i,(s,k) in enumerate(zip(S,'ABC')):
    a=fig.add_axes([0.03+i*0.155,0.04,0.145,0.27]); a.imshow(np.sqrt(s/s.max()),cmap=blue,interpolation='nearest')
    a.set_xticks([]);a.set_yticks([])
    for sp in a.spines.values(): sp.set_color('#1f3a66')
    a.set_title(f'slice_{k}  ({int((s>0).sum())} lit cells, max {s.max():g})',color='#9fbde8',fontsize=9)
# --- radial profiles
a=fig.add_axes([0.53,0.07,0.44,0.24],facecolor=bg)
for r,k,c in zip(R,'ABC',['#ff6b6b','#5fe08a','#5fa8ff']): a.plot(np.arange(64),r,color=c,lw=1.4,label='radial_'+k)
a.set_ylim(0,0.16); a.text(1.6,0.155,'r=1 peak 0.75 (all three)',color='#cfe3ff',fontsize=8,va='top')
a.legend(frameon=False,labelcolor='#cfe3ff',fontsize=8); a.tick_params(colors='#7f98bf',labelsize=8)
for sp in a.spines.values(): sp.set_color('#1f3a66')
a.set_xlabel('radius (cells from centre)',color='#9fbde8',fontsize=9)
fig.text(0.97,0.005,'50,752 B · sha256 3dd92a42…96c3aae · GGUF v3, 6 F32 tensors, read byte-exact',color='#5b7299',fontsize=8,ha='right')
fig.savefig(sys.argv[2],dpi=130,facecolor=bg)
