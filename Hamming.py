import numpy as np

def Hamming_distance_matrix(self, msa:np.array) -> np.array:
        "Fast implementation of generating a Hamming distance matrix"
        return (msa[:, None, :] != msa).sum(2)