import numpy as np 
from matplotlib import pyplot as plt 
from scipy.optimize import bisect
from scipy.io import FortranFile as FF

def get_thr(thr, densi, Npart):
    return (densi/thr).astype(int).sum() - Npart

#matplotlib settings
plt.rcParams["font.family"] = "Manjari"
#bold
plt.rcParams["font.weight"] = "bold"

#config ########################################
L = 100.0
Npart = 1e6
nx = 2**(int(round(np.log2(Npart**(1/3))))-1)
dx = L/nx
xvec = np.linspace(0+dx/2, L-dx/2, nx)
print(nx)

spectral_idx = -2.0
#######################################

#generate inhomogeneous distribution
freqs = np.fft.fftfreq(nx, d=dx)

kx, ky, kz = np.meshgrid(freqs, freqs, freqs)
k = np.sqrt(kx**2 + ky**2 + kz**2)

densi_tilde = np.random.normal(size=(nx, nx, nx)) + 1j * np.random.normal(size=(nx, nx, nx))
densi_tilde *= k**spectral_idx
densi_tilde[0, 0, 0] = 0.0
densi = np.fft.ifftn(densi_tilde).real

densi -= densi.min()
densi /= densi.sum()
densi *= Npart

thr = bisect(get_thr, densi.min()+1e-6, densi.max(), args=(densi, Npart))
parts = (densi/thr).astype(int)
print(f'Number of particles: {parts.sum()}')

xpart = [] 
ypart = []
zpart = []

for i in range(nx):
    for j in range(nx):
        for k in range(nx):
            npart = parts[i, j, k]
            if npart > 0:
                xpart += [np.random.uniform(xvec[i]-dx/2, xvec[i]+dx/2, npart)]
                ypart += [np.random.uniform(xvec[j]-dx/2, xvec[j]+dx/2, npart)]
                zpart += [np.random.uniform(xvec[k]-dx/2, xvec[k]+dx/2, npart)]

xpart = np.concatenate(xpart)
ypart = np.concatenate(ypart)
zpart = np.concatenate(zpart)

with FF('points.dat', 'w') as f:
    f.write_record(np.int32(len(xpart)))
    f.write_record(xpart.astype(np.float32))
    f.write_record(ypart.astype(np.float32))
    f.write_record(zpart.astype(np.float32)) 
###############################################


#plot#######################################
fig, ax = plt.subplots(figsize=(5,5), dpi = 150,)
ax.imshow(densi[:, :, int(nx/2)].T, cmap='PuBu_r', interpolation = 'gaussian', aspect='auto', origin='lower',)
#hide ticks
ax.set_xticks([])
ax.set_yticks([])
#Colorbar attached to the bottom
cbar = fig.colorbar(ax.images[0], ax=ax, orientation='horizontal', pad=0.002, )
cbar.set_label('Number of particles per pixel')
cbar.ax.xaxis.set_label_position('bottom')
plt.savefig('inhomo_dist.pdf', dpi=200, bbox_inches='tight', pad_inches=0.1)
#######################################

