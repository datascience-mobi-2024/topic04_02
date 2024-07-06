pos_corr = {'YR': 0.492442, 
                'RP': 0.480594, 
                'RG': 0.442533, 
                'R': 0.432905, 
                'WR': 0.392523, 
                'YP': 0.390709, 
                'LR': 0.386596, 
                'FR': 0.376553, 
                'VR': 0.370241, 
                'ER': 0.369174, 
                'RC': 0.357052, 
                'Rhelix': 0.349204, # percentage of R in all helices
                'MR': 0.343496, 
                'P': 0.340082, 
                'PG': 0.338434, 
                'EAmotif': 0.332202, # EA hinterenander
                'LP': 0.330337, 
                'EP': 0.328569, 
                'RH': 0.326193, 
                'AR': 0.324355, 
                'ARmotif': 0.3241, # AR hinterenander
                'NR': 0.323378}
neg_corr = {'QT': -0.52935, 
                'MQ': -0.507329, 
                'QS': -0.502697, 
                'QC': -0.493738, 
                'Q': -0.469765, 
                'QD': -0.466556, 
                'QH': -0.455041, 
                'NQ': -0.435562, 
                'IQ': -0.429683, 
                'FQ': -0.42363, 
                'WQ': -0.420057, 
                'QK': -0.41872, 
                'PolarAA': -0.406035, 
                'ST': -0.396185, 
                'Qhelix': -0.37895, 
                'TC': -0.364763, 
                'MT': -0.359458, 
                'TH': -0.346921, 
                'TD': -0.34469, 
                'SH': -0.327684, 
                'SC': -0.321686, 
                'T': -0.321208}

ideal_pos_value = {'AR': 0.19577028919092312, 
                       'VR': 0.16570774850840858,
                       'LR': 0.22057443760531548,
                       'LP': 0.19913083779456836,
                       'MR': 0.10189569954318338,
                       'FR': 0.1214807452476434,
                       'WR': 0.09732041140609597,
                       'NR': 0.10350733491352844,
                       'YR': 0.11505796446768121,
                       'YP': 0.09361436465693411,
                       'ER': 0.1750096588731577,
                       'EP': 0.15356605906241058,
                       'RH': 0.10505828276892132,
                       'RC': 0.0895071066024632,
                       'RP': 0.14867680267495723,
                       'RG': 0.1743851008086472,
                       'PG': 0.15294150099790005,
                       'PolarAA': 0.14402202391640712,
                       'R': 0.08506020124285214,
                       'P': 0.06361660143210504,
                       'ARmotif': 0.010758762261349177,
                       'EAmotif': 0.01633251571740483,
                       'Rhelix': 0.10373774208196636
                       }
ideal_neg_value = {'IQ': 0.05490211183163266,
                    'MQ': 0.03973022298244349,
                    'MT': 0.054568274035922376,
                    'FQ': 0.05931526868690351,
                    'WQ': 0.03515493484535609,
                    'NQ': 0.041341858352788544,
                    'QS': 0.057844351285310555,
                    'QT': 0.06062750041770345,
                    'QD': 0.06216042713743332,
                    'QH': 0.04289280620818142,
                    'QK': 0.06281534035966867,
                    'QC': 0.027341630041723304,
                    'ST': 0.07268240233878943,
                    'SH': 0.054947708129267414,
                    'SC': 0.03939653196280929,
                    'TD': 0.0769984781909122,
                    'TH': 0.0577308572616603,
                    'TC': 0.04217968109520218,
                    'Q': 0.022894724682112268,
                    'T': 0.03773277573559115,
                    'PolarAA': 0.14402202391640712, 
                    'Qhelix': 0.02831461743013483
                    }
#possible substitutions for each aminoacid, taken from literature
#https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003674
conserv_subst = {
        'A': ['D', 'E', 'G', 'S', 'T'],
        'C': ['G', 'R', 'S', 'W', 'Y'],
        'D': ['A', 'E', 'G', 'H', 'N', 'V', 'Y'],
        'E': ['A', 'D', 'G', 'K', 'Q', 'V'],
        'F': ['I', 'L', 'Y'],
        'G': ['A', 'C', 'D', 'E', 'R'],
        'H': ['D', 'L', 'N', 'P', 'Q', 'R', 'Y'],
        'I': ['F', 'L', 'M', 'N', 'V'],
        'K': ['E', 'M', 'N', 'Q', 'R', 'T'],
        'L': ['F', 'H', 'I', 'M', 'P', 'Q', 'R', 'V', 'W'],
        'M': ['I', 'K', 'L', 'R', 'T', 'V'],
        'N': ['D', 'H', 'I', 'K', 'S', 'T', 'Y'],
        'P': ['H', 'L', 'Q', 'R', 'S'],
        'Q': ['E', 'H', 'K', 'L', 'P', 'R'],
        'R': ['C', 'G', 'H', 'K', 'L', 'M', 'P', 'Q', 'T', 'W'],
        'S': ['A', 'C', 'N', 'P', 'T', 'W', 'Y'],
        'T': ['A', 'K', 'M', 'N', 'R', 'S'],
        'V': ['D', 'E', 'I', 'L', 'M'],
        'W': ['C', 'L', 'R', 'S'],
        'Y': ['C', 'D', 'F', 'H', 'N'],
        }
non_conservative_substitutions = {
        'A': ['P', 'V'],
        'C': ['F'],
        'F': ['C', 'S', 'V'],
        'G': ['S', 'V', 'W'],
        'I': ['K', 'R', 'S', 'T'],
        'K': ['I'],
        'L': ['S'],
        'P': ['A', 'T'],
        'Q': ['E', 'K'],
        'R': ['I', 'S'],
        'S': ['F', 'G', 'I', 'L', 'R'],
        'T': ['I', 'P'],
        'V': ['A', 'F', 'G'],
        'W': ['G'],
        }

AA_polar_neutral:list = ['N', 'Q', 'S', 'T', 'Y']

def rel_aa_comp(Sequence:str, AA_property): 
    count = 0
    for n in AA_property:
        for i in Sequence:  
            if n ==i:
                count += 1
    return count/len(Sequence)

#Clusters atoms if they are within a set distance
def cluster_calc(array, by_atom=False):
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
                atom_number.append(adjacency_matrix[list_components[0], n + 1])
                Cluster[f'Cluster {str(i)}'] = atom_number
    if by_atom:
        clust_inv = {}
        for cluster_name, atom_list in Cluster.items():
            for atom_number in atom_list:
                clust_inv[atom_number] = cluster_name
        return clust_inv
    else:
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
    import warnings
    
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    
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
    from ThERMOS import distance
    import numpy as np
    import warnings
    
    warnings.filterwarnings("ignore", category = RuntimeWarning)
    
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
    import numpy as np
    rows_with_nan = np.insert(np.array([np.all(np.isnan(array[1:, 1:]), axis=1)]),0, None) #find rows with all nan values
    cols_with_nan = np.insert(np.array([np.all(np.isnan(array[1:, 1:]), axis=0)]),0, None) #find columns with all nan values
    array = array[~rows_with_nan, :] #delete rows with all nan values
    array = array[:, ~cols_with_nan] #delete columns with all nan values
    array[:,0] = array[:,0].astype('int')
    return array

def pdb2AA(path, file_name, list_output=True):
    """
    Extracts the amino acid sequence from a PDB file.
    
    Args:
        path (str): Path to the directory containing the PDB file.
        file_name (str): Filename of the PDB file.

    Returns:
        list: A list containing the 1-letter amino acid sequence from the PDB file.
    """
    import os
    AA_dict = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
        "SEC": "U",
        "PYL": "O",
    }
    
    #extract aminoacid list
    with open(os.path.join(path, file_name), 'r') as f:
        ## path.join(path, file_name)
        import re
        aa1 = []
        aa_count = 1
        for line in f:
            line = line.replace('-', '  -')
            line = re.sub(r'([A])(\d)', r'\1 \2', line)
            if line.startswith('ATOM'):
                number = int(line.split()[5])
                aa = str(line.split()[3])
                if number == 1 and AA_dict[aa] not in aa1:
                    aa1.append(AA_dict[aa])
                if number == aa_count and len(aa1)+1 == number:
                    aa1.append(AA_dict[aa])
                aa_count = number
    if list_output:
        return aa1
    else:
        AA_string = ''.join(aa1)
        return AA_string

def AAA2A(AA:list, str=True):
    AA_dict = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
        "SEC": "U",
        "PYL": "O",
        }
    for n in range(len(AA)):
        AA[n] = AA_dict[AA[n]]
    if str:
        return "".join(AA)
    else:
        return AA
    
def ArraySlice (Array, possible_mutations):
    import numpy as np
    free_AA = Array
    index_list = [i.split('-')[1] for i in possible_mutations]
    mask = ~np.isin(Array, index_list)
    mask = np.all(mask, axis=1)
    filtered_AA = free_AA[mask]
    return filtered_AA

def fasparse(faspath):
    import numpy as np
    helixl = []
    sheetl = []
    data = np.loadtxt(faspath, dtype={'names': ('index', 'col1', 'col2', 'val1', 'val2', 'val3'),'formats': ('i4', 'S1', 'S1', 'f4', 'f4', 'f4')})
    indiceshelix = data['index'][data['col2'] == b'H']
    indicessheet = data['index'][data['col2'] == b'E']
    helix = np.split(indiceshelix, np.where(np.diff(indiceshelix) != 1)[0]+1)
    sheet = np.split(indicessheet, np.where(np.diff(indicessheet) != 1)[0]+1)
    return [helix, sheet]

def diff_weighted(feature_pos, feature_neg, aa:str, ideal_pos:dict, ideal_neg:dict, sec_prediction, sort = True, sum_only = False):
    """
    Calculates the weighted sum of deviations for a given amino acid ('aa') based on positive and negative features.

    Args:
        feature_pos (list): List representing the positive features.
        feature_neg (List): List representing the negative features.
        aa (str): String of the amino acid sequence.
        ideal_pos (dict): Dictionary containing ideal values for positive features.
        ideal_neg (dict): Dictionary containing ideal values for negative features.
        sort (bool, optional): Flag indicating whether to sort the features by deviation (default: True).

    Returns:
        tuple: A tuple containing two elements:
            - sum_dev (float): The total weighted sum of deviations for all features.
            - sorted_keys (list, optional): If sort=True, a list of features sorted by their weighted deviation (highest first).
                - WT_weight (dict, optional): If sort=False, a dictionary containing the weighted deviation for each feature.
    """
    from ThERMOS import rel_aa_comp
    import operator

    AA_polar = 'NQSTY'
    
    WT_weight = {}
    sum_dev = 0
    
    #positive features
    for key in feature_pos: 
        if len(key) <=2:
            
            weighted_diff = abs(rel_aa_comp(aa, key) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
            
        elif 'motif' in key:
            weighted_diff = abs(rel_aa_comp(aa, key[0:2]) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
    
        elif 'helix' in key:
            feature = key[0]
            helices = sec_prediction[0]
            aa_helix = []
            for aa_pos in range(len(aa) + 1):
                for helix in helices:
                    if aa_pos+1 in helix:
                        aa_helix.append(aa[aa_pos])
            if len(aa_helix) > 0:
                weighted_diff = abs(rel_aa_comp(''.join(aa_helix), feature) - (ideal_pos[key] * feature_pos[key]))
                WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
    
    #negative features
    for key in feature_neg:
        if len(key) <=2:
            weighted_diff = abs(rel_aa_comp(aa, key) - ideal_neg[key] * feature_neg[key])
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff
            
        elif 'motif' in key:
            weighted_diff = abs(rel_aa_comp(aa, key[0:2]) - (ideal_pos[key] * feature_pos[key]))
            WT_weight[key] = weighted_diff
            sum_dev += weighted_diff    
            
        elif 'Polar' in key:
            weighted_diff = abs(rel_aa_comp(aa, AA_polar) - (ideal_neg[key] * feature_neg[key]))
        
        elif 'helix' in key:
            feature = key[0]
            helices = sec_prediction[0]
            aa_helix = []
            for aa_pos in range(len(aa) + 1):
                for helix in helices:
                    if aa_pos+1 in helix:
                        aa_helix.append(aa[aa_pos])
            if len(aa_helix) > 0:
                weighted_diff = abs(rel_aa_comp(''.join(aa_helix), feature) - (ideal_neg[key] * feature_neg[key]))
                WT_weight[key] = weighted_diff
            sum_dev += weighted_diff

    if sum_only == True:
        return sum_dev
    elif sort:
        sorted_keys = sorted(WT_weight.items(), key=operator.itemgetter(1), reverse=True)
        return sum_dev, sorted_keys
    else:
        return sum_dev, WT_weight #sum_dev is a positive value of all deviations, higher values indicate a worse fit

def mut_apply(AA_list, Mut_list):
    """
    Applies a list of mutations to a list of amino acids.

    Args:
        AA_list (list): A list containing the original amino acid sequence.
        Mut_list (list): A list of mutation strings defining the substitutions to be applied in the form ('WT-POS-MUT')

    Returns:
        list: A new list containing the amino acid sequence after applying the mutations.
    """
    if len(Mut_list) > 0:
        if 'M' in Mut_list[0]:
            Mut_list = Mut_list[1:]
        for n in Mut_list:
            AA_pos = int(n.split('-')[1]) - 1
            AA_mut = n.split('-')[2]
            AA_list[AA_pos] = AA_mut
        return AA_list
    else:
        return(print('Mut list is empty'))

def mut_live_test (AA_list, Mut_list, pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction):
    from ThERMOS import mut_apply
    from ThERMOS import diff_weighted
    AAs = ''.join(AA_list)
    WT_sum, WT_diff = diff_weighted(pos_corr, neg_corr, AAs, ideal_pos_value, ideal_neg_value, sec_prediction)
    AA_mut = mut_apply(AA_list, Mut_list)
    AA_muts = ''.join(AA_mut)
    MUT_sum, WT_diff = diff_weighted(pos_corr, neg_corr, AA_muts, ideal_pos_value, ideal_neg_value, sec_prediction)
    Diff = abs(WT_sum) - abs(MUT_sum) # calculates the difference between the weighted sum of deviations before and after the mutation
    
    # if the difference is positive the mutation is beneficial
    # if Diff is negative the mutation is not beneficial
    return Diff #The higher the difference the better the mutation
    
def mutator_rand(AAs_list, substitutions, threshhold = 100, seed = 0):
    from itertools import product
    import random
    count = 0
    random.seed(seed)
    
    keys = list(substitutions.keys())
    random.shuffle(keys)
    # Create a list of lists, each containing tuples of (position, substitution)
    while count < threshhold:

        subst_options = [[(pos, subst) for subst in [AAs_list[int(pos)-1]] + substitutions[pos]] for pos in keys]
        if not subst_options:
            break
        
        for combination in product(*subst_options):
            # Start with the original protein sequence
            prot_variation = list(AAs_list)
            # Apply each substitution in the combination
            for pos, subst in combination:
                prot_variation[int(pos)-1] = subst
            # Yield the new protein variation as a string
            yield ''.join(prot_variation)
            count += 1
            
            if count >= threshhold:
                break
                    
def functional_aa(input_path, pdb_file, output_path, df=False):
    """
    This function selects atoms involved in various interactions (salt bridges, hydrogen bonds, 
    and van der Waals interactions) from a protein structure.

    Args:
        input_path (str): Path to the directory containing the protein structure file.
        pdb_file (str): Filename of the protein structure file in PDB format.
        output_path (str): Path to the directory where the output PQR file will be saved.
        df (bool, optional): If True, returns a pandas DataFrame containing the information. 
                                Defaults to False (returns a NumPy array).

    Returns:
        np.ndarray | pd.DataFrame: A NumPy array containing the selected atom information 
                                   or a pandas DataFrame if `df` is True.
    """
    
    #import necessary functions
    from ThERMOS import salt_bridge
    from ThERMOS import H_bond_calc
    from ThERMOS import VdW_interaction
    from ThERMOS import pdb2pqr
    from ThERMOS import remove_nan
    import os
    import numpy as np
    import pandas as pd
    import re
    
    # get protein name and pqr file name
    prot_name = pdb_file.split('-')[1]
    
    #create pqr file
    pqr_file = f'{(pdb_file.split('.')[0]).split('-')[1]}.pqr'

    if os.path.isfile(os.path.join(output_path, f'{(pdb_file.split(".")[0]).split("-")[1]}.pqr')):
        print('Pqr file already exists')
    else:    
        pdb2pqr(input_path, output_path, pdb_file)

    # Calculate atom features
    Salt_bridge = salt_bridge(input_path, pdb_file)
    H_bond = H_bond_calc(output_path, pqr_file)
    VdW_clust, VdW_vol = VdW_interaction(input_path, pdb_file, by_atom = True)
    
    # extract the values for the proteins from the dictionary and delete atoms that dont have a feature (if applicable)
    Salt_bridge = remove_nan(Salt_bridge[prot_name])
    H_bond = remove_nan(H_bond[prot_name][:,:,0])

    VdW_clust = VdW_clust[prot_name]

    #create lists with all aminoacid that are part of a feature
    atom_S =list(Salt_bridge[0,1:])
    atom_HA = list(H_bond[0,1:])
    atom_HD = list(H_bond[1:,0])

    # creates an atom_dict that contains the atom number and the feature it is part of
    atom_dict = {}
    for lst, identifier in [(atom_S, "Salt_bridge"), (atom_HA, "Hbond_acc"), (atom_HD, "Hbond_don")]:
        for atom_number in lst:
            if atom_number in atom_dict:
                atom_dict[atom_number] = [atom_dict[atom_number], identifier]
            else:
                atom_dict[atom_number] = identifier
    # Add van der Waals interaction information to the dictionary
    for k,v in VdW_clust.items():
        if k in atom_dict:
            atom_dict[k] = [atom_dict[k], v]
        else: atom_dict[k] = v
    atom_sorted = {k: atom_dict[k] for k in sorted(atom_dict)}

    # create a dataframe with the atom number and the feature it is part of
    prot_df = pd.DataFrame(columns = ['Protein','Aminoacid','Aminoacid_number', 'Atom_number', 'Feature'])
    Protein_array = np.empty((0, 5))
    with open (os.path.join(output_path, pqr_file)) as f:
        prot_df_list = []
        for line in f:
            #line = line.replace('-', '  -')
            #line = re.sub(r'([A])(\d)', r'\1 \2', line)
            if line.startswith('ATOM'):
                atom_number = int(line.split()[1])
                feature = atom_sorted.get(atom_number)
                if atom_number in atom_sorted.keys():
                    atom_line = np.array([[str(prot_name),str(line.split()[3]), int(line.split()[4]),int(line.split()[1]), str(feature)]])
                    Protein_array = np.append(Protein_array, atom_line, axis=0)
                    
    # return the dataframe if df is True
    if df:
        Prot_df = pd.DataFrame(Protein_array, columns = ['Protein','Aminoacid','Aminoacid_number', 'Atom_number', 'Feature'])
        return Prot_df
    else:
        return Protein_array

def SASA_calc(path, pdb_files=None):
    from Bio.PDB import PDBParser
    from Bio.PDB.SASA import ShrakeRupley
    import os
    SASA_dict = {}
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files]
    for pdb_file in pdb_files:
        struct = PDBParser(QUIET=1).get_structure(pdb_file.split('-')[1], os.path.join(path, str(pdb_file)))
        ShrakeRupley().compute(struct, level = 'S')
        SASA_dict[pdb_file.split('-')[1]] = struct.sasa
    return SASA_dict

def pdb2pqr(input_path, output_path,pdb_files=None):
#INFO:Please cite:  Jurrus E, et al.  Improvements to the APBS biomolecular solvation software suite.  Protein Sci 27 112-128 (2018).
#INFO:Please cite:  Dolinsky TJ, et al.  PDB2PQR: expanding and upgrading automated preparation of biomolecular structures for molecular simulations. Nucleic Acids Res 35 W522-W525 (2007).
    import os
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(input_path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files]
    for pdb_file in pdb_files:
        name = f'{(pdb_file.split('.')[0]).split('-')[1]}.pqr'
        os.system(f'pdb2pqr "{os.path.join(input_path, str(pdb_file))}" "{os.path.join(output_path, name)}" -ff={'AMBER'} --noop')

def salt_bridge(path, pdb_files=None): #https://www.bioinformation.net/003/002800032008.pdf
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    from ThERMOS import distance
    import re
    
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files]
    Salt_bridges = dict()
    
    for pdb_file in pdb_files:
        Asp_Glu_array = np.empty((0, 4))
        Lys_Arg_His_array = np.empty((0, 4))
        
        with open(os.path.join(path, str(pdb_file))) as f:
            for line in f:
                #line = line.replace('-', '  -')
               #line = re.sub(r'([A])(\d)', r'\1 \2', line)
                if line.startswith('ATOM'):
                    if ('ASP' in line and 'OD' in line) or ('GLU' in line and 'OE' in line):
                        line_array = np.array([[line[7:12].strip(), line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]])
                        line_array = line_array.astype('float64')
                        Asp_Glu_array = np.append(Asp_Glu_array, line_array, axis = 0)
                    if ('LYS' in line and 'NZ' in line) or ('ARG' in line and 'NH' in line) or ('HIS' in line and 'NE' in line) or ('HIS' in line and 'ND' in line):
                        line_array = np.array([[line[7:12].strip(), line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]])
                        line_array = line_array.astype('float64')
                        Lys_Arg_His_array = np.append(Lys_Arg_His_array, line_array, axis = 0)

            Salt_bridges[str(pdb_file).split('-')[1]] = distance(Asp_Glu_array, Lys_Arg_His_array, 4)
    return Salt_bridges

def VdW_interaction(path, pdb_files=None, by_atom = False):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    import re
    
    from ThERMOS import distance
    from ThERMOS import cluster_calc
    from ThERMOS import intersect_vol
    
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files]
    VdW_cluster = {}
    VdW_volume = {}

    for pdb_file in pdb_files:
        with open(os.path.join(path, str(pdb_file))) as f:
            Atom_array = np.empty((0, 4))
            #VdW_radii = {'C': 3.1, 'N': 2.95, 'O': 2.96} # Van der Waals radii in Angstrom enlarged by watermolecule radius 1.4 A (https://academic.oup.com/nar/article/49/W1/W559/6279848#267025710)
            #C_C = 6.2
            #C_N = 6.05
            #C_O = 6.06
            #N_N = 5.9
            #N_O = 5.91
            #O_O = 5.92
            X_array = np.array([[]])
            Atom_list = ['C', 'N', 'O']
            for line in f:
                #line = line.replace('-', '  -')
                #line = re.sub(r'([A])(\d)', r'\1 \2', line)
                line = line.strip()
                if line.startswith('ATOM'):
                    if ('LEU' in line and line[12:17].strip() in Atom_list) or ('VAL' in line and line[12:17].strip() in Atom_list) or ('ILE' in line and line[12:17].strip() in Atom_list):
                        line_array = np.array([[line[7:12].strip(), line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]])
                        line_array = line_array.astype('float64')
                        Atom_array = np.append(Atom_array, line_array, axis = 0)
                        if 'C' in line[12:17]:
                            X_array = np.append(X_array, int('0'))
                        elif 'N' in line[12:17]:
                            X_array = np.append(X_array, int('1'))
                        elif 'O' in line[12:17]:
                            X_array = np.append(X_array, int('2'))
            
            Atom_distance = distance(Atom_array, Atom_array, 6, remove_nan = False)
            
            Atom_distance = np.nan_to_num(Atom_distance)
            VdW_cluster[str(pdb_file).split('-')[1]] = cluster_calc(Atom_distance, by_atom)
            
            Atom_distance_nan = np.where(Atom_distance==0, np.nan, Atom_distance)
            Atom_volume = intersect_vol(Atom_distance_nan, 6, 6)
            VdW_volume[str(pdb_file).split('-')[1]] = Atom_volume
            
    return VdW_cluster, VdW_volume
                
def H_bond_calc(path, pqr_files=None):
    #https://www.sciencedirect.com/science/article/pii/S2665928X20300246?via%3Dihub#sec2
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    from ThERMOS import distance
    from ThERMOS import angle_calc
    
    Donor_dict = {'GLN': [('NE2', 'HE21'), ('NE2', 'HE22')], 
                'GLU': [('OE2', 'HE2')], 
                'ASP': [('OD2', 'HD2')], 
                'ASN': [('ND2', 'HD21'), ('ND2', 'HD22')], 
                'HIS': [('NE2', 'HE2'), ('ND1', 'HD2')], 
                'LYS': [('NZ', 'HZ1'), ('NZ', 'HZ2'), ('NZ', 'HZ3')], 
                'ARG': [('NE', 'HE'), ('NH1', 'HH11'), ('NH1', 'HH12'), ('NH2', 'HH21'), ('NH2', 'HH22')], 
                'SER': [('OG', 'HG')], 
                'THR': [('OG1', 'HG1')], 
                'TRP': [('NE1', 'HE1')], 
                'TYR': [('OH', 'HH')]}

    Acceptor_dict = {'GLN': [('OE1')], 
                'ASP': [('OD1'), ('OD2')], 
                'ASN': ['OD1'], 
                'GLU': [('OE1'), ('OE2')], 
                'SER': ['OG'], 
                'THR': ['OG1']}
    HB_dict = {}
    if pqr_files is None:
        pqr_files = [f for f in os.listdir(path) if f.endswith('.pqr')]
    if isinstance(pqr_files, str):
        pqr_files = [pqr_files]
    for pqr_file in pqr_files:
        with open(os.path.join(path, str(pqr_file))) as f:
            Donor_array = np.empty((0, 4))
            H_array = np.empty((0, 4))
            Acceptor_array = np.empty((0, 4))
            aa_cache = []
            atom_cache = []
            for line in f:
                line = line.replace('-', '  -')
                if line.startswith('ATOM'):
                    if not aa_cache:
                        aa_cache.append(line.split()[3])
                        aa_cache.append(line.split()[4])
                        atom_cache.append(line)
                    elif aa_cache[1] == line.split()[4]:
                        atom_cache.append(line)
                    elif aa_cache[1] != line.split()[4]:
                        if aa_cache[0] in Donor_dict.keys():
                            sub_donor = Donor_dict[aa_cache[0]] #extracts list (with tupels of Donor, Hydrogen) for the amino acid
                            for n in sub_donor:
                                donor_match = [entry for entry in atom_cache if n[0] in entry]
                                h_match = [entry for entry in atom_cache if n[1] in entry]
                                if donor_match and h_match:
                                    d_line = np.array([[int(donor_match[0].split()[1]), float(donor_match[0].split()[5]), float(donor_match[0].split()[6]), float(donor_match[0].split()[7])]]) 
                                    h_line = np.array([[int(h_match[0].split()[1]), float(h_match[0].split()[5]), float(h_match[0].split()[6]), float(h_match[0].split()[7])]])
                                    Donor_array = np.append(Donor_array, d_line, axis=0)
                                    H_array = np.append(H_array, h_line, axis=0)
                        if aa_cache[0] in Acceptor_dict.keys():
                            sub_acc = Acceptor_dict[aa_cache[0]] #extracts list of acceptors for the amino acid
                            for n in sub_acc:
                                acc_match = [entry for entry in atom_cache if n in entry]               
                                a_line = np.array([[int(acc_match[0].split()[1]), float(acc_match[0].split()[5]), float(acc_match[0].split()[6]), float(acc_match[0].split()[7])]])
                                Acceptor_array = np.append(Acceptor_array, a_line, axis=0)
                        aa_cache = []
                        atom_cache = [] 

        angle = angle_calc(Donor_array, H_array, Acceptor_array)
        HB_dict[str(pqr_file).split('.')[0]] = angle

    return HB_dict
                                            
def AA2s4pred (directory_S4pred, output_path, AA_seq, prot, remove_file = None):
    import os
    from ThERMOS import fasparse
    os.getcwd()
    # call s4pred and create fas file
    fastapath = os.path.join(output_path,f'{prot}.fasta')
    faspath = os.path.join(output_path,f'{prot}.fas')
    abs_fasta = os.path.abspath(fastapath)
    abs_fas = os.path.abspath(faspath)

    if os.path.isfile(fastapath):
        print(f'fasta file already exists')
    else:
        with open(os.path.join(output_path,f'{prot}.fasta'), "w") as fasta_file:
            fasta_file.write(f">{prot}\n{AA_seq}\n")
        
    if os.path.isfile(faspath):
        print('fas file already exists')
    else:
        os.chdir(directory_S4pred)
        os.system(f'python3 run_model.py "{abs_fasta}" > "{abs_fas}"')
        os.chdir('../../')
     
    #read fas file and   
    sec_pred = fasparse(abs_fas)
    if remove_file:
        os.remove(fastapath)
        os.remove(faspath) 

    return sec_pred 

def free_aa (path, pdb_file, functional_aa):
    """
    Identifies and collects free amino acids from a PDB file.

    This function takes a path to a PQR file, the filename of the PQR file, and a NumPy array containing protein information as input.
    It iterates through the PQR file and identifies residues that are not involved in Salt bridges, Hydrogen bonds, and Van der Waals interactions.

    Args:
        path (str): Path to the directory containing the PQR file.
        pqr_file (str): Filename of the PQR file.
        prot_arr (np.ndarray): NumPy array containing protein information (assumed to have residue types in the 2nd column).

    Returns:
        np.ndarray: A NumPy array containing information about free amino acids (protein name, residue name, residue number).
    """
       
    import os
    import numpy as np
    import re
    AA_dict = {
        "ALA": "A",
        "CYS": "C",
        "ASP": "D",
        "GLU": "E",
        "PHE": "F",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LYS": "K",
        "LEU": "L",
        "MET": "M",
        "ASN": "N",
        "PRO": "P",
        "GLN": "Q",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "VAL": "V",
        "TRP": "W",
        "TYR": "Y",
        "SEC": "U",
        "PYL": "O",
    }
    functional_aa = sorted(set(functional_aa[:,2]))
    free_aa = np.empty((0, 3))
    prot_name = pdb_file.split('-')[1]
    with open(os.path.join(path, str(pdb_file))) as f:
        for line in f:
            line = line.replace('-', '  -')
            line = re.sub(r'([A])(\d)', r'\1 \2', line)
            if line.startswith('ATOM'):
                aa_number = line.split()[5]
                if aa_number not in functional_aa and aa_number not in free_aa[:,2]:
                    aa_line = np.array([[str(prot_name), AA_dict[line.split()[3]], line.split()[5]]])
                    free_aa = np.append(free_aa, aa_line, axis=0)
    return free_aa                

def mutator_rational(AA_list:list, free_AA, deviation, pos_corr:dict, neg_corr:dict, conserv_substitution, ideal_pos_value, ideal_neg_value, cutoff, sec_prediction):
    """
    Generates a list of potential mutations based on deviations and correlations.

    Args:
        AAs (str): The amino acid sequence.
        free_AA (np.ndarray): Array containing information about free amino acids Col1: Prot name, Col2: aminoacid position, Col3: Aminoacid (one letter code).
        deviation (list): List containing deviations from ideal values for features.
        pos_corr_list (list): List of amino acids that positively correlate with desired features.
        sorted_freq_pos (list): Possibly sorted list of frequencies for amino acids contributing to positive features (usage unclear).
        neg_corr_list (list): List of amino acids that negatively correlate with desired features.
        conserv_substitution (dict): Dictionary containing a list of possible conservative substitutions for each amino acid.

    Returns:
        list: Mutated protein as a list
        list: List of mutations (AA-POS-AA)
    """
    
    import itertools
    from ThERMOS import rel_aa_comp
    from ThERMOS import mut_apply
    from ThERMOS import mut_live_test
    from ThERMOS import AA_polar_neutral
    
    pos_corr_list = list(pos_corr.keys())
    neg_corr_list = list(neg_corr.keys())
    AA_mut_list = AA_list
    AAs = ''.join(AA_list)
    mut_AAs = AAs
    free_AA_dict = {a: b for a, b in zip(free_AA[:,2], free_AA[:,1] )} # create dictionary from array the key is the absolute aminoacid position and value is the aminoacid
    mut_list = []
    first_entry = deviation[0][0]
    max_increase = 1.3 # maximum increase of relative amino acid composition
    Diff_start = cutoff
    helices = sec_prediction[0]
    
    

    # determine possible substitutions if the first entry is a single amino acid
    if len(first_entry) == 1:
    
        if first_entry in pos_corr_list: #checks if amount of aminoacid should be increaed
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                AA_subst = conserv_substitution[key]
                
                
                for k in pos_corr_list: # Check if the current amino acid is in the positive correlation list
                    if len(k) == 1:     #selects the first feature with one aminoacid
                        if rel_aa_comp(mut_AAs, k) < ideal_pos_value[k]:    #checks if the relative composition is suboptimal
                            if k in AA_subst:
                                mut_aa = (aminoacid + '-' + key + '-' + k)
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                                if Diff > Diff_start:
                                    mut_list.append(mut_aa)
                                    AA_mut_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                    mut_AAs = ''.join(AA_mut_list)  #needs to be adjustet to not change original sequence
                                    break


        elif first_entry in neg_corr_list: #checks if amount of aminoacid should be decreased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                if aminoacid == first_entry:
                    AA_subst = conserv_substitution[key]
                    
                    #tries to substitute the aminoacid to the aminoacids that comes first in the sorted_freq_pos list
                    for k in pos_corr_list:
                        if aminoacid == k:
                            if rel_aa_comp(mut_AAs, k) < ideal_pos_value[k]:
                                mut_aa = aminoacid + '-' + key + '-' + k
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                                if Diff > Diff_start: # live test if mutation is benefitical
                                    AA_mut_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                    mut_list.append(mut_aa)
                                    mut_AAs = ''.join(AA_mut_list)
                                    break
                        
                        

    # determine possible substitutions if the first entry is a pair of amino acids
    elif len(first_entry) == 2:
        
        if first_entry in pos_corr_list: # checks if the amount of aminoacids should be increased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                #if aminoacid not in sorted_freq_pos: #checks if the aminoacid overall positively contributes to one of the pos_corr features
                AA_subst = conserv_substitution[key] # list of possible substitutions for the current amino acid
                
                #mutation
                for n in range(len(deviation)): # iterates through all deviations and takes the first deviation (from pos corr) which can be increaed
                    entry = deviation[n][0]
                    if entry in pos_corr_list:
                        subst = [] #creates a list of possible substitutions that increase one of the amino acids in the highest entry that positively correlates
                        for k in AA_subst:
                            if k in first_entry:
                                subst = k
                                
                        if len(subst) == 1: # if only one substitution increases one of the aminoacids this substitution will be used
                            mut_aa = aminoacid + '-' + key + '-' + subst
                            Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                            mut_rel = rel_aa_comp(AAs, subst)
                            mut_rel_max = rel_aa_comp(AAs, subst) * max_increase
                            if Diff > Diff_start and mut_rel <= mut_rel_max:
                                AA_mut_list = mut_apply(AA_list, [mut_aa]) #applies mutation to the protein sequence
                                mut_list.append(mut_aa)
                                mut_AAs = ''.join(AA_mut_list)
                                break
                            
                        elif len(subst) == 2: # if the aminoacid can be substituted to both aminoacids in first_entry choose the one that has the least frequency
                            comp = [(aa, rel_aa_comp(AAs, aa)) for aa in subst] # calculate the relative amino acid composition of the possible substitutions
                            lowest_comp = min(comp, key=lambda pair: pair[1]) # find the amino acid with the lowest relative composition
                            mut_aa = aminoacid + '-' + key + '-' + lowest_comp[0]
                            mut_rel = rel_aa_comp(AAs, lowest_comp[0])
                            mut_rel_max = rel_aa_comp(AAs, lowest_comp[0]) * max_increase
                            Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                            if Diff > Diff_start and mut_rel <= mut_rel_max:
                                AA_mut_list = mut_apply(AA_list, [mut_aa])
                                mut_list.append(mut_aa) # append the mutation to the mutation list
                                mut_AAs = ''.join(AA_mut_list)
                                break

        
                    
        elif first_entry in neg_corr_list: # checks if the amount of aminoacids should be decreased
            for key in free_AA_dict:
                aminoacid = free_AA_dict[key]
                if aminoacid in first_entry: #checks if the aminoacid is present in the first_entry
                    AA_subst = conserv_substitution[key] # list of possible substitutions for the current amino acid
                    
                    #mutation               
                    for n in range(len(deviation)): # iterates through all deviations and takes the first deviation (from pos corr) which can be increaed
                        entry = deviation[n][0]
                        if entry in pos_corr_list: 
                            subst = [] #creates a list of possible substitutions that increase one of the amino acids in the first_entry
                            for k in AA_subst:
                                if k in entry:
                                    subst.append(k)
                            if len(subst) == 1: # if only one substitution increases one of the aminoacids this substitution will be used
                                mut_aa = aminoacid + '-' + key + '-' + subst[0]
                                mut_rel = rel_aa_comp(mut_AAs, subst[0])
                                mut_rel_max = rel_aa_comp(AAs, subst[0]) * max_increase
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                                if Diff > Diff_start and mut_rel <= mut_rel_max:
                                    AA_mut_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa)
                                    mut_AAs = ''.join(AA_mut_list) 
                                    break
                                
                            elif len(subst) == 2: # if the aminoacid can be substituted to both aminoacids in first_entry choose the one that has the least frequency
                                comp = [(aa, rel_aa_comp(mut_AAs, aa)) for aa in subst] # calculate the relative amino acid composition of the possible substitutions
                                lowest_comp = min(comp, key=lambda pair: pair[1]) # find the amino acid with the lowest relative composition
                                mut_aa = aminoacid + '-' + key + '-' + lowest_comp[0]
                                
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                                mut_rel = rel_aa_comp(mut_AAs, lowest_comp[0])
                                mut_rel_max = rel_aa_comp(AAs, lowest_comp[0]) * max_increase
                                
                                if Diff > Diff_start and mut_rel <= mut_rel_max:
                                    AA_mut_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa) # append the mutation to the mutation list
                                    mut_AAs = ''.join(AA_mut_list)
                                    break
    

    elif 'motif' in first_entry:
        if first_entry in pos_corr_list: # checks if the amount of aminoacids should be increased, as of know motifs are only in the pos_corr_list
            for key in free_AA_dict:
                prev_key = int(key) -2
                curr_key = int(key) -1
                for_key = int(key) 
                
                
                #handles edge-case for first aminoacid in the sequence
                if key == 1: #if the aminoacid is the first in the sequence
                    curr_aa = AAs[curr_key]
                    for_aa = AAs[for_key]
                                                                    
                    if curr_aa not in first_entry[0] and for_aa in first_entry[1]:
                        curr_subst = conserv_substitution[key]
                        for s in curr_subst:
                            if s in first_entry[0]:
                                s_rel = rel_aa_comp(mut_AAs, s)
                                s_rel_max = rel_aa_comp(AAs, s) * max_increase
                                if s in first_entry[0] and s_rel <= s_rel_max:
                                    mut_aa = f'{curr_key}-{for_key}-{s}'
                                    Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                                    if Diff > Diff_start:
                                        AA_mut_list = mut_apply(AA_list, [mut_aa])
                                        mut_list.append(mut_aa)
                                        mut_AAs = ''.join(AA_mut_list)
                                        break
                
                
                #handles edge-case if for the last aminoacid in the sequence
                elif for_key >= len(AAs):  #if the aminoacid is the last in the sequence
                    prev_aa = AAs[prev_key]
                    curr_aa = AAs[curr_key]
                    if prev_aa not in first_entry[0] or curr_aa not in first_entry[1]:
                        
                        if prev_aa in first_entry[0] and curr_aa not in first_entry[1]:
                            curr_subst = conserv_substitution[key]
                            for s in curr_subst:
                                if s in first_entry[1]:
                                    s_rel = rel_aa_comp(mut_AAs, s)
                                    s_rel_max = rel_aa_comp(AAs, s) * max_increase
                                    if s in first_entry[1] and s_rel <= s_rel_max:
                                        mut_aa = f'{prev_key}-{curr_key}-{s}'
                                        Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                                        if Diff > Diff_start:
                                            AA_mut_list = mut_apply(AA_list, [mut_aa])
                                            mut_list.append(mut_aa)
                                            mut_AAs = ''.join(AA_mut_list)
                                            break
                
                
                else:            
                    aa_back = AAs[prev_key]    # aminoacid before the current aminoacid
                    aa_current = AAs[curr_key]   # current aminoacid
                    aa_for = AAs[for_key]     # aminoacid after the current aminoacid
                    
                    if aa_back in first_entry[0] and aa_current not in first_entry[1]: #checks if previous aminoacid is part of the motif, and current aminoacid not
                        curr_subst = conserv_substitution[key] # possible substitutions for the current aminoacid
                        for s in curr_subst:
                            s_rel = rel_aa_comp(mut_AAs, s)
                            s_rel_max = rel_aa_comp(AAs, s) * max_increase
                            if s in first_entry[1] and s_rel <= s_rel_max:
                                mut_aa = f'{curr_key}-{for_key}-{s}'
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction) # live test if mutation is benefitical
                                if Diff > Diff_start:
                                    AA_mut_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa)
                                    mut_AAs = ''.join(AA_mut_list)
                                    break
                    
                    if aa_current not in first_entry[0] and aa_for in first_entry[1]: #checks if current aminoacid is not part of the motif, and next aminoacid is
                        curr_subst = conserv_substitution[key] # possible substitutions for the current aminoacid
                        for s in curr_subst:
                            s_rel = rel_aa_comp(mut_AAs, s)
                            s_rel_max = rel_aa_comp(AAs, s) * max_increase
                            if s in first_entry[0] and s_rel <= s_rel_max:
                                mut_aa = f'{curr_key}-{for_key}-{s}'
                                Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                                if Diff > Diff_start:
                                    AA_mut_list = mut_apply(AA_list, [mut_aa])
                                    mut_list.append(mut_aa)
                                    mut_AAs = ''.join(AA_mut_list)
                                    break
                                
    elif 'Polar' in first_entry:
        if first_entry in neg_corr_list:
            for key in free_AA_dict:
                        aminoacid = free_AA_dict[key]
                        if aminoacid in AA_polar_neutral:
                            aa_subst = conserv_substitution[key]
                            best_Diff = mut_live_test(AA_mut_list, [aminoacid], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                            for subst in aa_subst:
                                if subst not in AA_polar_neutral:
                                    mut_aa = f'{aminoacid} + - + {key} + - + {subst}'
                                    Diff = mut_live_test(AA_mut_list, [subst], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                                    if Diff > best_Diff:
                                        best_Diff = Diff
                                        best_mut_aa = subst
                            AA_mut_list = mut_apply(AA_list, [best_mut_aa])
                            mut_list.append(best_mut_aa)
                            mut_AAs = ''.join(AA_mut_list)
                    
    elif 'helix' in first_entry:
        if first_entry in pos_corr_list:
            aa2increase = first_entry[0]         
            for pos in free_AA_dict:
                for helix in helices:
                        if first_entry in helix:
                            aminoacid = free_AA_dict[pos]
                            if aminoacid != aa2increase:
                                poss_subst = conserv_substitution[pos]
                                for subst in poss_subst:
                                    if subst == aa2increase:
                                        mut_aa = f'{aminoacid}-{pos}-{subst}'
                                        Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, sec_prediction)
                                        if Diff > Diff_start:
                                            AA_mut_list = mut_apply(AA_list, [mut_aa])
                                            mut_list.append(mut_aa)
                                            mut_AAs = ''.join(AA_mut_list)
                                            break
        
        if first_entry in neg_corr_list:
            aa2reduce = first_entry[0]
            for pos in free_AA_dict:
                for helix in helices:
                    if first_entry in helix:
                        aminoacid = free_AA_dict[pos]
                        if aminoacid == aa2reduce:
                            poss_subst = conserv_substitution[pos]
                            
                            best_Diff = mut_live_test(AA_mut_list, [aminoacid], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value)
                            for subst in poss_subst:
                                if subst != aa2reduce: #probably redundant, bcs if aminoacid is aa2reduce it is not in poss_subst
                                    mut_aa = f'{aminoacid}-{pos}-{subst}'
                                    Diff = mut_live_test(AA_mut_list, [mut_aa], pos_corr, neg_corr, ideal_pos_value, ideal_neg_value)
                                    if Diff > best_Diff:
                                        best_Diff = Diff
                                        best_subst = subst
                                        best_mut_aa = mut_aa
                            AA_mut_list = mut_apply(AA_list, [best_subst])
                            mut_list.append(best_mut_aa)
                            mut_AAs = ''.join(AA_mut_list)
                            
                            

                             
    return AA_mut_list, mut_list

def Subst_reducer(sec_pred:list, conserv_subst_dict:dict, free_AA_dict:dict, seed):
    """
    Reduces the possible substitutions for each amino acid based on secondary structure predictions.

    Args:
        sec_pred: A list of length 2 containing secondary structure predictions for each position in the sequence.
          - sec_pred[0]: List of characters representing helix ('H') or coil ('-') predictions for each position.
          - sec_pred[1]: List of characters representing sheet ('E') or coil ('-') predictions for each position.
        conserv_subst_dict: A dictionary where keys are amino acids and values are lists of their conservative substitutions.
        free_AA_dict: A dictionary where keys are amino acid position (in protein) and value is aminoacid.

    Returns:
        A dictionary where keys are amino acids and values are reduced lists of possible substitutions based on secondary structure predictions.
    """
    import random
    
    random.seed(seed)
    
    helix_forming = ['E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H', 'I', 'W', 'F']
    sheet_forming = ['L', 'M', 'V', 'I', 'Y', 'C', 'W', 'F', 'T', 'U']
    
    #helixvalues = {'E':1.59,'A':1.41,'L':1.34,'M':1.3,'Q':1.27,'K':1.23,'R':1.21,'H':1.05,'V':0.9,'I':1.09,'Y':0.74,'C':0.66,'W':1.02,'F':1.16,'T':0.76,'G':0.43,'N':0.76,'P':0.34,'S':0.57,'D':0.99,'U':0.66}
    #sheetvalues = {'E':0.52,'A':0.72,'L':1.22,'M':1.14,'Q':0.98,'K':0.69,'R':0.84,'H':0.8,'V':1.87,'I':1.67,'Y':1.45,'C':1.4,'W':1.35,'F':1.33,'T':1.17,'G':0.58,'N':0.48,'P':0.31,'S':0.96,'D':0.39,'U':1.4}


    helix = sec_pred[0]
    sheet = sec_pred[1]
    Possible_subst = {}
    
    for key in free_AA_dict:
        aminoacid = free_AA_dict[key]
        if any(key in n for n in helix):
            possible_subst= list(set(conserv_subst_dict[aminoacid]).intersection(set(helix_forming)))
        elif any(key in n for n in sheet):
            possible_subst = list(set(free_AA_dict[key]).intersection(set(sheet_forming)))
        else:
            possible_subst = conserv_subst_dict[aminoacid]
        
        possible_subst.append(aminoacid)
            
        random.shuffle(possible_subst)
        Possible_subst[key] = possible_subst
        
    return Possible_subst

def ThERMOS(pdb_path, pdb_file, pqr_output_path, locked_aa_pos=None, Deep_mut=True, iterations=100, cutoff_value = -0.005, threshhold = 10000, seed = 0, remove_files = None):
    """
    This function performs protein mutation analysis to improve the thermal stability of a protein, with minimal changes to the structure
    (as measured by melting point). It takes a PDB file path, filename, and path for PQR output as input.

    Args:
        pdb_path (str): Path where the pdb file is stored.
        pdb_file (str): Name of the PDB file.
        pqr_output_path (str): Path where pqr file will be saved (also fasta and fas).
        locked_aa_pos (list, optional): List of amino acid positions that should not be mutated (defaults to None).
        Deep_mut (bool, optional): Flag whether to include random mutation in screening
            (random + rational, defaults to True).
        iterations (int, optional): Number of iterations for rational improvement (defaults to 100).
        cutoff_value (float, optional): Cutoff value for mutation selection in rational improvement 
            (defaults to -0.005).
        threshhold (int, optional): Threshold for random mutation acceptance (higher value leads to more mutations, 
            defaults to 10000).
        seed (int, optional): Seed for the random number generator (ensures reproducibility, defaults to 0).

    Returns:
        list: A list containing three elements:
            - Tuple: (WT_SPARC object, best_SPARC object) - Wild-type and best mutated protein SPARC predictions.
            - Tuple: (WT amino acid list, best mutated protein amino acid list) - Amino acid sequences.
            - Tuple: (WT deviation sum, best mutated protein deviation sum) - Deviations of the protein structures.
    """
    #import functions

    from ThERMOS import AA2s4pred
    from ThERMOS import diff_weighted
    from ThERMOS import mutator_rand
    from ThERMOS import mutator_rational
    from ThERMOS import functional_aa
    from ThERMOS import free_aa
    from ThERMOS import Subst_reducer
    from ThERMOS import pdb2AA
    from ThERMOS import ArraySlice
    from ThERMOS import pos_corr, neg_corr, ideal_pos_value, ideal_neg_value, conserv_subst, non_conservative_substitutions

    from SPARC import SPARC
    
    from heapq import heappop, heappush
    import heapq    
    import os
   

    #extract protein features
    aa_list = pdb2AA(pdb_path, pdb_file)
    aa_locked = functional_aa(pdb_path, pdb_file, pqr_output_path)
    aa_free = free_aa(pdb_path, pdb_file, aa_locked)
    aa_str = ''.join(aa_list)
    free_AA_dict = {a: b for a, b in zip(aa_free[:,2], aa_free[:,1] )} # create dictionary from array the key is the absolute aminoacid position and value is the aminoacid
    
    if locked_aa_pos is not None:
        for pos in locked_aa_pos:
            if free_AA_dict.get(pos) is not None:
                free_AA_dict.pop(pos)
    
    sec_prediction = AA2s4pred('./data/s4pred', pqr_output_path, aa_str, pdb_file, remove_file = remove_files)
    possible_substitutions = Subst_reducer(sec_prediction, conserv_subst, free_AA_dict, seed = seed)

    #calculate WT deviations
    WT_dev_sum, WT_dev = diff_weighted(pos_corr, neg_corr, aa_list, ideal_pos_value, ideal_neg_value, sec_prediction)

    if Deep_mut:
        #randoly mutate protein (within given constraints), get top 10 mutations
        top_variations = []
        largest_variations = []
        heappush(top_variations, (WT_dev_sum, WT_dev, aa_str))  # Placeholder for lowest score
        heappush(largest_variations, (0, WT_dev, aa_str))  # Placeholder for largest variation
        
        Mut_seq_str = mutator_rand(aa_list, possible_substitutions, threshhold = threshhold, seed = seed)
        for Mut_prot in Mut_seq_str:
            Mut_dev_sum, Mut_dev = diff_weighted(pos_corr, neg_corr, Mut_prot, ideal_pos_value, ideal_neg_value, sec_prediction, sort = True)
            Str_dev = sum(c1 != c2 for c1, c2 in zip(''.join(Mut_prot), aa_str))

            #save the top 10 scores
            if len(top_variations) < 11:
                heappush(top_variations, (Mut_dev_sum, Mut_dev, Mut_prot)) 
                    
            elif Mut_dev_sum < top_variations[0][0]:
                heappush(top_variations, (Mut_dev_sum, Mut_dev, Mut_prot))
                heappop(top_variations)
            
            #saves top 10 largest variations
            if len(largest_variations) <11:
                heappush(largest_variations, (Str_dev, Mut_dev, Mut_prot))
            elif Str_dev > largest_variations[0][0]:
                heappush(largest_variations, (Str_dev, Mut_dev, Mut_prot))
                heappop(largest_variations)
                
        heappush(top_variations, (WT_dev_sum, WT_dev, aa_str))
        Random_creation = list(heapq.merge(top_variations, largest_variations))
        print('Random mutation finished')
        
        
    #define variables for iteration
    prev_Mut_prot_list = aa_list    
    prev_Mut_dev = WT_dev           
    aa_available = aa_free          
    
    #define variables for best protein#
    if Deep_mut:
        best_Mut_prot_list = list(heappop(top_variations)[2])
        best_Mut_dev_sum = heappop(top_variations)[0]
        best_Mut_dev = heappop(top_variations)[1]
        best_aa_available = aa_available
    else:
        best_Mut_prot_list = list(aa_list)
        best_Mut_dev_sum = WT_dev_sum
        best_Mut_dev = WT_dev
        best_aa_available = aa_available
        Random_creation = [(WT_dev_sum, WT_dev, aa_str)]
    
    # Initiate top 5 best variations of rational improvement to calculate melt point
    top_top_variations = []
    heappush(top_top_variations, (float(WT_dev_sum), WT_dev, aa_str))
    
    
    #use top 10 mutated sequences and use rational improvement      
    best_iteration = 0 #(used to track how many iterations are necessary, currently ~2-3 seems best)
    for mut_seq in Random_creation:
        prev_Mut_prot_list = list(mut_seq[2])
        prev_Mut_dev_sum = mut_seq[0]
        prev_Mut_dev = mut_seq[1]

        for k in range(iterations):
            Mut_prot_list, possible_mutations = mutator_rational(
                                                    AA_list = prev_Mut_prot_list, 
                                                    free_AA = aa_available, 
                                                    deviation = prev_Mut_dev,
                                                    pos_corr = pos_corr, 
                                                    neg_corr =  neg_corr, 
                                                    conserv_substitution = possible_substitutions,
                                                    ideal_pos_value = ideal_pos_value, 
                                                    ideal_neg_value = ideal_neg_value,
                                                    cutoff = cutoff_value,
                                                    sec_prediction = sec_prediction
                                                    ) #f_value = cutoff, calculates list of possible mutations
            
            Mut_dev_sum, Mut_dev = diff_weighted(pos_corr, neg_corr, Mut_prot_list, ideal_pos_value, ideal_neg_value, sec_prediction) # calculate deviation of mutated protein sequence
            
            best_dev = float(top_top_variations[0][0])
            if len(top_top_variations) < 6:
                heappush(top_top_variations, (float(Mut_dev_sum), Mut_dev, Mut_prot_list))

            elif float(Mut_dev_sum) < best_dev:
                heappush(top_top_variations, (Mut_dev_sum, Mut_dev, Mut_prot_list))
                heappop(top_top_variations)
                
            if abs(best_Mut_dev_sum) > abs(Mut_dev_sum):
                best_Mut_prot_list = Mut_prot_list
                best_Mut_dev_sum = Mut_dev_sum
                best_Mut_dev = Mut_dev  
                best_possible_mutations = possible_mutations #get list of best mutations (AA-POS-AA), depreciated, bcs random mutator doesn't output this
                best_iteration = str(k+1)
                #aa_available = ArraySlice(aa_available, possible_mutations) #updates available aminoacids, so that each aminoacid can only be mutated once
                
            elif Mut_dev[0][0] == prev_Mut_dev[0][0] and abs(Mut_dev[0][1]-prev_Mut_dev[0][1]) < 0.001:
                break
                
            #update variables for next iteration            
            prev_possible_mutations = possible_mutations    #list of mutations (AA-POS-AA) (prev_possible_mutations can be printed if needed)
            prev_Mut_prot_list = Mut_prot_list                        #Mutated protein as a list with one AA per entry
            prev_Mut_dev = Mut_dev
            prev_Mut_dev_sum = Mut_dev_sum
    
        
    #for top_hit in top_top_variations:
    #-------SPARC implementation missing--------#
    wt_sparc = SPARC(aa_str, pdb_file.split('-')[1], './data', './data/s4pred')
    best_temp = wt_sparc[0]
    wt_temp = wt_sparc[0]  
    best_SPARC = wt_sparc
    
    #select best mutation based on melt point
    for top in top_top_variations:
        top_SPARC = SPARC(''.join(top[2]), pdb_file.split('-')[1], './data', './data/s4pred')
        top_temp = top_SPARC[0]
        if top_SPARC[0] > best_temp:
            best_temp = top_SPARC[0]
            best_SPARC = top_SPARC
            best_Mut_prot_list = top[2]
            best_Mut_dev_sum = top[0]
            best_Mut_dev = top[1]
    
    
    
    Improvement = WT_dev_sum - best_Mut_dev_sum
    
    return [(wt_sparc, best_SPARC), (aa_list, best_Mut_prot_list), (WT_dev_sum, best_Mut_dev_sum)]

def ThERMless(mut_temp, wt_temp, wt_protein, mut_protein, name, cutoff = 0.9, min_diff = 0, sec_prediction=None, fast=False, ):
    from ThERMOS import diff_weighted
    from ThERMOS import pos_corr, neg_corr, ideal_pos_value, ideal_neg_value
    from ThERMOS import AA2s4pred
    from SPARC import SPARC
    from heapq import heappop, heappush, nlargest, heapify

    import os
    #diff_weighted(feature_pos, feature_neg, aa:str, ideal_pos:dict, ideal_neg:dict, sec_prediction, sort = True, sum_only = False):
    #define variables for iteration
    start_Tm_diff = mut_temp-wt_temp
    wt_str = ''.join(wt_protein)
    best_mut_str = ''.join(mut_protein)
    Tm_diff = start_Tm_diff

            
    if sec_prediction == None:
        sec_prediction = AA2s4pred('./data/s4pred', './data', wt_str, name)
        os.remove(os.path.join('./data', f'{name}.fasta'))
        os.remove(os.path.join('./data', f'{name}.fas'))
    
    wt_diff = diff_weighted(pos_corr, neg_corr, wt_str, ideal_pos_value, ideal_neg_value, sec_prediction, sum_only = True)
    start_mut_diff = diff_weighted(pos_corr, neg_corr, best_mut_str, ideal_pos_value, ideal_neg_value, sec_prediction, sum_only = True)
    mut_diff = start_mut_diff

    count  =  0
    while (Tm_diff >= start_Tm_diff * cutoff) and (Tm_diff >= min_diff):
        sorted_wt_screen = [(float('inf'), 'dummy')]
        for i in range(len(best_mut_str)):
            if best_mut_str[i] != wt_str[i]:
                mut_str = best_mut_str[:i] + wt_str[i] + best_mut_str[i+1:]
                mut_diff = diff_weighted(pos_corr, neg_corr, mut_str, ideal_pos_value, ideal_neg_value, sec_prediction, sum_only = True) #calculate deviation to fully mutated protein
                heappush(sorted_wt_screen, (mut_diff, mut_str))
        count  +=1

        #get new mutated protein with least difference to original mutated protein 
        best_mut_diff = heappop(sorted_wt_screen)[0]
        best_mut_str = heappop(sorted_wt_screen)[1]
        sparc_screen = SPARC(best_mut_str, name, './data', './data/s4pred')
        Tm_diff = sparc_screen[0][0] - wt_temp
        No_mutations = 0 
        for i in range(len(best_mut_str)):
            if best_mut_str[i] != wt_str[i]:
                No_mutations += 1
        
        if best_mut_str == wt_str:
            break
    
    return (best_mut_str, best_mut_diff, sparc_screen[0][0])