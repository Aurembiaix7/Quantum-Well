import numpy as np
meff=0.06 # Massa efectiva de la banda de conducció [units of electron rest mass] (Material Parameters)
hbar2om=7.62 # hbar^2/electron rest mass [eV Amg^2] (Physical Parameters)
VB=2.3 # Alçada de la barrera [eV]

# Simulation parameters
nombrePunts=600
print(nombrePunts)
delta=1. # Separació entre punts [Amg]
width=60. # Well Width [Amg]

# Setting up Hamiltonian Matrix
H=np.zeros((nombrePunts,nombrePunts))
offdiag=-hbar2om/(2.,meff*delta**2)
for i in range (nombrePunts): # Offdiagonal Hamiltonian elements
    H[i,i+1]=offdiag
    H[i,i-1]=offdiag
    H[i,i]=-2.*offdiag
    nwidth=width/delta # Quants punts tindrà la barrera
    half=nombrePunts//2 # Divisió entera per nombrePunts imparell
    if i<half-nwidth//2 or i>half+nwidth//2:
        H[i,i]+=VB # Afegir la Valence Band

# Diagonalise Hamiltioian Matrix
energies, eigenstates=np.linalg.eigh(H)

# Plot
import matplotlib.pyplot as plt
plt.figure()
x=np.arange(nombrePunts)*delta # Per tenir unitats de distància
plt.plot(x,eigenstates[:,0])

np.set_printoptions(precision=3)
print(energies)

print(eigenstates)

print("Programa acabat")