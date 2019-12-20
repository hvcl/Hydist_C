import numpy as np
NANGDAY = 0.0

h = np.loadtxt("bandodosau.txt")
h = np.transpose(h)
N, M = h.shape
h = np.flip(h, 1)
h = np.pad(h, ((1, 2), (1, 2)), 'edge')
h = h + NANGDAY

np.savetxt("bandodosau.txt",h)


hsnham = np.loadtxt("hsnham.txt")
hsnham = np.transpose(hsnham)
hsnham = np.flip(hsnham, 1)
hsnham = np.pad(hsnham, ((1, 2), (1, 2)), 'edge')
np.savetxt("hsnham.txt",hsnham)
# VISCOINDX: hsnhot
VISCOINDX = np.loadtxt("hsnhotroiA.txt")
VISCOINDX = np.transpose(VISCOINDX)
VISCOINDX = np.flip(VISCOINDX, 1)
VISCOINDX = np.pad(VISCOINDX, ((1, 2), (1, 2)), 'edge')
np.savetxt("hsnhotroiA.txt",VISCOINDX)

# hs ma sat day
Fw = np.loadtxt("Fw_map.txt")
Fw = np.transpose(Fw)
Fw = np.flip(Fw, 1)
Fw = np.pad(Fw, ((1, 2), (1, 2)), 'edge')
np.savetxt("Fw_map.txt",Fw)

