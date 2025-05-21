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



def display_hamiltonian(matrix):
    labels = [r'$xy \uparrow$', r'$xy \downarrow$', r'$yx \uparrow$', r'$yx \downarrow$', 
              r'$z^2 \uparrow$', r'$z^2 \downarrow$', r'$xz \uparrow$', r'$xz \downarrow$', 
              r'$(x^2 - y^2) \uparrow$', r'$(x^2 - y^2) \downarrow$']

    mat_color = 'cornflowerblue'
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot each element in the matrix with a black border and display the values
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1, edgecolor='black', 
                                 linewidth=1, facecolor=mat_color)
            ax.add_patch(rect)
            if matrix[i][j] == 0:
                ax.text(j, i, '0', ha='center', va='center', fontsize=12)
            else:
                ax.text(j, i, f'{matrix[i][j]:.3f}', ha='center', va='center', fontsize=12)

    # Set the limits and aspect ratio of the plot
    ax.set_xlim(-0.5, 9.5)
    ax.set_ylim(9.5, -0.5)
    ax.set_aspect('equal')

    # Set ticks
    ax.set_xticks(np.arange(10))
    ax.set_yticks(np.arange(10))
    ax.set_xticklabels(labels, fontsize=12, rotation=45)
    ax.set_yticklabels(labels, fontsize=12)
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False)


    # Display the plot
    plt.show()

hamiltonian = convert_to_nparr("H_10x10.dat")
print("Full 10x10 Hamiltonian")
print(hamiltonian)

# display_hamiltonian(hamiltonian)

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

# # Labels for d-orbitals
# d_orbitals = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2y2']
# 
# contributions = np.abs(sorted_eigenvectors)**2
# normalized_contributions = contributions / np.sum(contributions, axis=0)
# 
# fig, ax = plt.subplots(figsize=(10, 6))
# 
# x = np.arange(1, len(sorted_eigenvalues) + 1)
# 
# # Plot each orbital's contribution as a stacked bar
# bottom = np.zeros(len(sorted_eigenvalues))  # To stack the bars
# for i in range(len(d_orbitals)):
#     ax.bar(x, normalized_contributions[i, :], bottom=bottom, label=d_orbitals[i])
#     bottom += normalized_contributions[i, :]
# 
# # Set xtick labels to be the eigenvalues rounded to 4 decimal places
# eigenvalue_labels = [f"{val:.4f}" for val in sorted_eigenvalues]
# ax.set_xticks(x)  # Keep the ticks evenly spaced
# ax.set_xticklabels(eigenvalue_labels, fontsize=12)
# 
# # Set plot labels and title
# ax.set_xlabel('Eigen-energy', fontsize=20)
# ax.set_ylabel('Normalized Contribution', fontsize=20)
# ax.set_title('Normalized Contribution of Each d-Orbital to Eigenvectors', fontsize=14)
# ax.set_xticks(x)
# ax.set_ylim(0, 1)
# ax.legend(title='d-Orbitals', fontsize = 15)
# 
# # Display the plot
# plt.tight_layout()
# plt.show()

##################

# Labels for d-orbitals
d_orbitals = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2y2']

contributions = np.abs(sorted_eigenvectors)**2
normalized_contributions = contributions / np.sum(contributions, axis=0)

tableau_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

# Left subplot: Stacked bar plot for normalized contributions
x = np.arange(1, len(sorted_eigenvalues) + 1)
bottom = np.zeros(len(sorted_eigenvalues))  # To stack the bars
for i in range(len(d_orbitals)):
    ax1.bar(x, normalized_contributions[i, :], bottom=bottom, label=d_orbitals[i], color=tableau_colors[i])
    bottom += normalized_contributions[i, :]

eigenvalue_labels = [f"{val:.4f}" for val in sorted_eigenvalues]
ax1.set_xticks(x)
ax1.set_xticklabels(eigenvalue_labels, fontsize=12)

# Set plot labels and title for the left subplot
ax1.set_xlabel('Eigen-energy', fontsize=20)
ax1.set_ylabel('Normalized Contribution', fontsize=20)
ax1.set_title('Normalized Contribution of Each d-Orbital to Eigenvectors', fontsize=14)
ax1.set_ylim(0, 1)
ax1.legend(title='d-Orbitals', fontsize=15)

# Right subplot: Plot eigenvalues with colored line segments for each orbital
y_positions = np.arange(len(sorted_eigenvalues))  # Generate positions for each eigenvalue
y_positions = sorted_eigenvalues  # Generate positions for each eigenvalue

# Loop through each eigenvalue
for j, eigenvalue in enumerate(sorted_eigenvalues):
    # Start at 0.25 and increment based on contributions
    x_start = 0.25
    x_end = 0.75
    total_length = x_end - x_start

    # Loop through each orbital's contribution to the current eigenvalue
    for i in range(len(d_orbitals)):
        contribution_length = total_length * normalized_contributions[i, j]  # Proportional length
        ax2.hlines(y=y_positions[j], xmin=x_start, xmax=x_start + contribution_length, 
                   colors=tableau_colors[i], lw=4)
        x_start += contribution_length  # Move the starting point for the next segment

    ax2.text(0.76, y_positions[j], f"{eigenvalue:.4f}", va='center', fontsize=12, color='black')

ax2.set_title('Eigenvalues with Orbital Contributions', fontsize=14)
ax2.set_xlim(0, 1)
ax2.set_xticks([])
ax2.grid(True)

plt.tight_layout()
plt.show()
