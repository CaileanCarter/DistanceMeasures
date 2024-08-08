import numpy as np
import blosum as bl

def BLOSUM_distance(arr1, arr2, BLOSUM) -> float:
        "Calculate the BLOSUM distance between two sequences"
        return sum((BLOSUM[a1][a2] for a1,a2 in zip(arr1,arr2)))

def BLOSUM_distance_matrix( kmer_groups, shape):
        blos_dist = np.full(shape, np.nan) # can change default value
        blos_dist_func = self.BLOSUM_distance
        bl_matrix = self.bl_matrix
        indeces = {k : i for i,k in enumerate(self.kmers)} 

        for cl in kmer_groups:
            if len(cl) == 1:
                continue
            for arr1, arr2 in combinations(cl, 2):
                arr1_i = indeces[arr1]
                arr2_i = indeces[arr2]
                blos_dist[arr1_i, arr2_i] = blos_dist[arr2_i, arr1_i] = blos_dist_func(arr1, arr2, bl_matrix)
        return blos_dist

def BLOSUM_rescale(arr):
        "Rescale BLOSUM distance to a range of 0-1"
        shape = arr.shape
        arr_flat = arr.flatten() 
        notnas = arr_flat[~np.isnan(arr_flat)]
        _max = max(notnas)
        _min = min(notnas)

        arr_flat = np.nan_to_num(arr_flat, nan=_min)
        arr_flat = 1-((arr_flat-_min)/(_max-_min))
        return arr_flat.reshape(shape)