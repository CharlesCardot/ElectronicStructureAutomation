import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

def convert_to_nparr(filename):
    mat = open(Path.cwd() / filename).readlines()
    mat = [x.replace("{","") for x in mat[0].split("}")][0:-2]
    for key,val in enumerate(mat):
        mat[key] = val[1:] if val[0] == "," else val

    mat = np.asarray([[float(y) for y in x.split(",")] for x in mat])
    return mat

hamiltonian = convert_to_nparr("H_10x10.dat")
print("Full 10x10 Hamiltonian")
print(hamiltonian)

mat = hamiltonian[::2, ::2] # Remove spin
print("\n5x5 Hamiltonian, no spin")
print(mat)

# Calculate eigenvalues and eigenvectors
eigenvalues, eigenvectors = np.linalg.eigh(mat)

# eigenvectors = eigenvectors.T

# Sort eigenvalues and corresponding eigenvectors (lowest to highest energy)
sorted_indices = np.argsort(eigenvalues)
sorted_eigenvalues = eigenvalues[sorted_indices]
sorted_eigenvectors = eigenvectors[:, sorted_indices]

# Print sorted eigenvalues and eigenvectors
print("\nEigenvalues (from lowest to highest energy):")
print(sorted_eigenvalues)

print("\nCorresponding eigenvectors:")
print(sorted_eigenvectors)  

#############
# Making Plot
#############

# Labels for d-orbitals
d_orbitals = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2y2']

contributions = np.abs(sorted_eigenvectors)**2
normalized_contributions = contributions / np.sum(contributions, axis=0)

fig, ax = plt.subplots(figsize=(10, 6))

x = np.arange(1, len(sorted_eigenvalues) + 1)

# Plot each orbital's contribution as a stacked bar
bottom = np.zeros(len(sorted_eigenvalues))  # To stack the bars
for i in range(len(d_orbitals)):
    ax.bar(x, normalized_contributions[i, :], bottom=bottom, label=d_orbitals[i])
    bottom += normalized_contributions[i, :]

# Set xtick labels to be the eigenvalues rounded to 4 decimal places
eigenvalue_labels = [f"{val:.4f}" for val in sorted_eigenvalues]
ax.set_xticks(x)  # Keep the ticks evenly spaced
ax.set_xticklabels(eigenvalue_labels, fontsize=12)

# Set plot labels and title
ax.set_xlabel('Eigen-energy', fontsize=20)
ax.set_ylabel('Normalized Contribution', fontsize=20)
ax.set_title('Normalized Contribution of Each d-Orbital to Eigenvectors', fontsize=14)
ax.set_xticks(x)
ax.set_ylim(0, 1)
ax.legend(title='d-Orbitals')

# Display the plot
plt.tight_layout()
plt.show()
