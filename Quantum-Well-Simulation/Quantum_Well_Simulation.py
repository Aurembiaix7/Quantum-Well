import numphy as np
meff=0.06 # Massa efectiva de la banda de conducció [units of electron rest mass] (Material Parameters)
hbar2om=7.62 # hbar^2/electron rest mass [eV Amg^2] (Physical Parameters)
VB=2.3 # Alçada de la barrera [eV]

# Simulation parameters
N # nombre de punts del grid
N=600
delta=1. # Separació entre punts [Amg]
width=60. # Well Width [Amg]

# Setting up Hamiltonian Matrix
H=np.zero((N, N))
offdiag=-hbar2om(2*delta**2*meff)
for i in range (N-2): # Offdiagonal Hamiltonian elements
H[i,i+1]=offdiag
H[i,i-1]=offdiag
for i in range(N-1): # Anirà de 0 a N-1
H[i,i]=-2.*offdiag
# Afegim el terme de la barrera
nwidth=width/delta # Quants punts tindrà la barrera
half=N//2 # Divisió entera per N imparell
for i in range(N-1): # Anirà de 0 a N-1
if (i<half-nwidth//2)or(i>half+nwidth//2)
H[i,i]+=VB # Afegir la Valence Band

# Diagonalise Hamiltioian Matrix
energies, eigenstates=np.linalg.eigh(H)

# Plot
import matplotlib.pyplot as plot
plt.figure()
x=np.arange(N)*delta # Per tenir unitats de distància
plt.plot(x,eigenstates[:,0])
