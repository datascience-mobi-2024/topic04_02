########################### Used for 3D structure analysis ###########################

#Clusters atoms if they are within a set distance
def cluster_calc(array):
    Cluster ={}
    from scipy.sparse.csgraph import connected_components
    from scipy.sparse import lil_matrix
    adjacency_matrix = lil_matrix((int(len(array)), int(len(array))))
    for i in range(len(array)):
        adjacency_matrix[0,i] = array[0,i]
        adjacency_matrix[i,0] = array[i,0]
    for i in range(1, len(array)):
        for j in range(i+1, len(array)):
            overlap_volume = array[i, j]
            adjacency_matrix[i, j] = overlap_volume

    adjacency_matrix = adjacency_matrix.tocsr()
    labels, n_components = connected_components(adjacency_matrix[1:,1:])

    set_components = set(n_components) #unique number of clusters
    list_components = n_components.tolist()

    for i in range(len(set_components)):
        atom_number = []
        for n in range(len(list_components)):
            if list_components[n] == i:
                atom_number.append(adjacency_matrix[list_components[0],n+1])   
                Cluster[f'Cluster {str(i)}'] = atom_number
    return Cluster

#Calculates the volume of intersection between two spheres
def intersect_vol(array, r1:str, r2:str):
    import numpy as np
    r1 = float(r1)
    r2 = float(r2)
    vol = array
    d = array[1:,1:]
    vol[1:,1:] = (np.pi*(r2+r1-d)**2 * (d**2 + 2*d*r1 - 3*r1**2 + 2*d*r2 + 6*r1*r2 - 3*r2**2))/(12*d)
    return vol

#Calculates distance based on two arrays, with set cutoff as a maximum distance allowed
def distance (array1, array2, cutoff, remove_nan=True):
    from scipy.spatial.distance import cdist
    import numpy as np
    distance = cdist(array1[:,1:], array2[:,1:], metric='euclidean') #calculate distance
    distance = np.concatenate((np.array([array2[:,0]]), distance), axis=0) #add atom number from  array2
    distance = np.concatenate((np.insert(np.array([array1[:,0]]), 0, None).reshape(-1,1), distance), axis=1) #add atom number from array1
    distance[1:, 1:][distance[1:, 1:] > cutoff] = np.nan #set distance > cutoff to nan
    if remove_nan == False:
        return distance
    elif remove_nan == True:
        rows_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=1)]),0, None)
        rows_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=1)]),0, None) #find rows with all nan values
        cols_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=0)]),0, None) #find columns with all nan values
        distance = distance[~rows_with_nan, :] #delete rows with all nan values
        distance = distance[:, ~cols_with_nan] #delete columns with all nan values
        return distance
    else:
        raise ValueError('remove_nan must be either True or False')