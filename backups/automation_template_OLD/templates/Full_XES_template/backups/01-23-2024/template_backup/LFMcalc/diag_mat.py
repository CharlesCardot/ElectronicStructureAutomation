import numpy as np
from scipy.linalg import orth

from pathlib import Path

def convert_to_nparr(filename):
    mat = open(Path.cwd() / filename).readlines()
    mat = [x.replace("{","") for x in mat[0].split("}")][0:-2]
    for key,val in enumerate(mat):
        mat[key] = val[1:] if val[0] == "," else val

    mat = np.asarray([[float(y) for y in x.split(",")] for x in mat])
    return mat

mat = convert_to_nparr("H_10x10.dat")
print(mat)

e_values = np.round(np.linalg.eig(mat)[0],8)
e_vectors = np.round(np.linalg.eig(mat)[1].T,8)
print("Eigenvlaues: ", e_values)

pairs = [np.where(e_values == e_val)[0].tolist() for e_val in e_values]
degen_e_vals = []
_ = [degen_e_vals.append(x) for x in pairs if x not in degen_e_vals]
print("Degenerate Eigenvalues: ", degen_e_vals)


# Check for degenerate subspaces
print(e_vectors)
for sublist in degen_e_vals:
    if not((len(sublist) % 2) == 0):
        print("ERROR: Somehow you have an odd number of states in your degenerate subspace")
    else:
        print("Found spatial degeneracy in e_vectors, finding orthonormal basis within degenerate subspace...")
        temp = np.asarray([e_vectors[x] for x in sublist])
        orth_temp = np.round(orth(temp.T),8)
        for key,x in enumerate(sublist):
            e_vectors[x] = orth_temp.T[key]
        

# Write e_vectors to file,
# keeping degenerate eigenvectors grouped
newline_count = 0
with open("new_opps.txt", "w") as f:
    for key,group in enumerate(degen_e_vals):
        for val in group:
            print(val)
            print(e_values[val])
            print(", ".join([str(x) for x in e_vectors[val]]))
            f.write(", ".join([str(x) for x in e_vectors[val]]))
            newline_count += 1
            if newline_count < len(e_vectors):
                f.write("\n")





