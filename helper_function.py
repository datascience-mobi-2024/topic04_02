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
def intersect_vol(array, r1:str, r2:str, dis = None):   
    import numpy as np
    r1 = float(r1)
    r2 = float(r2)
    vol = array
    if dis is not None:
        d = dis
    else:
        d = array[1:,1:]
    vol[1:,1:] = (np.pi*(r2+r1-d)**2 * (d**2 + 2*d*r1 - 3*r1**2 + 2*d*r2 + 6*r1*r2 - 3*r2**2))/(12*d)
    return vol

#Calculates distance based on two arrays, with set cutoff as a maximum distance allowed
def distance (array1, array2, cutoff = None, remove_nan=True):
    from scipy.spatial.distance import cdist
    import numpy as np
    distance = cdist(array1[:,1:], array2[:,1:], metric='euclidean') #calculate distance
    distance = np.concatenate((np.array([array1[:,0]]).T, distance), axis=1) #add atom number from  array1
    distance = np.concatenate(((np.insert(np.array([array2[:,0]]), 0, None).reshape(-1,1)).T, distance), axis=0) #add atom number from array1
    if cutoff is not None:
        distance[1:, 1:][distance[1:, 1:] >= cutoff] = np.nan #set distance > cutoff to nan
        
    if remove_nan == False:
        return distance
    elif remove_nan == True:
        rows_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=1)]),0, None) #find rows with all nan values
        cols_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=0)]),0, None) #find columns with all nan values
        distance = distance[~rows_with_nan, :] #delete rows with all nan values
        distance = distance[:, ~cols_with_nan] #delete columns with all nan values
        distance[:,0] = distance[:,0].astype('int')
        return distance
    else:
        raise ValueError('remove_nan must be either True or False')
    
    
def angle_calc(Donor_array, H_array, Acceptor_array): #https://www.sciencedirect.com/science/article/pii/S2665928X20300246?via%3Dihub
    from helper_function import distance
    import numpy as np
    d_DH = np.full((0,2,2), fill_value = np.nan)
    for n in range(len(Donor_array)):
        DH_temp = distance(np.array([Donor_array[n,:]]), np.array([H_array[n,:]]))
        d_DH = np.concatenate((d_DH, DH_temp.reshape((1,) + DH_temp.shape)), axis=0)
    d_HA = distance(H_array, Acceptor_array, cutoff = 3.5,remove_nan=False)
    d_DA = distance(Donor_array, Acceptor_array, remove_nan=False)
    angle = np.full((d_DA.shape[0], d_DA.shape[1]), fill_value = np.nan)
    angle[0,:] = d_DA[0,:].T
    for n in range(d_DH.shape[0]):
            theta = np.arccos((d_DH[n,1,1]**2 + d_HA[n,1:]**2 - d_DA[1:,1:][n,:]**2)/(2*d_DH[n,1,1]*d_HA[1:,1:][n,:]))
            theta[(theta < 100* np.pi/180) | (theta > np.pi)] = np.nan
            angle[n,0] = d_DA[n,0]
            angle[1:,1:][n,:] = theta 
    H_id = np.full((angle.shape[0], angle.shape[1]), fill_value = np.nan)
    H_id[1:,0] = d_DH[:,1,0]
    angle = np.dstack((angle, H_id))
    return angle



def remove_nan(array):
    from numpy import np
    rows_with_nan = np.insert(np.array([np.all(np.isnan(array[1:, 1:]), axis=1)]),0, None) #find rows with all nan values
    cols_with_nan = np.insert(np.array([np.all(np.isnan(array[1:, 1:]), axis=0)]),0, None) #find columns with all nan values
    array = array[~rows_with_nan, :] #delete rows with all nan values
    array = array[:, ~cols_with_nan] #delete columns with all nan values
    array[:,0] = array[:,0].astype('int')
