#function that calculates the relative amino acid composition based on sequence entry and definded amino acid property
def rel_aa_comp(Sequence:str, AA_property): 
    count = 0
    for n in AA_property:
        for i in Sequence:  
            if n ==i:
                count += 1
    return count/len(Sequence)

def rel_aa(Sequence:str, AA_property:str) -> str: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count

########################### Used for 3D structure analysis ###########################
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

#https://www.bioinformation.net/003/002800032008.pdf
def salt_bridge(path, pdb_files=None):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
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

            from helper_function import distance
            Salt_bridges[str(pdb_file).split('-')[1]] = distance(Asp_Glu_array, Lys_Arg_His_array, 4)
    return Salt_bridges

def VdW_interaction(path, pdb_files=None, by_atom = False):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    import re
    
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
            
            from helper_function import distance
            Atom_distance = distance(Atom_array, Atom_array, 6, remove_nan = False)
            
            from helper_function import cluster_calc
            Atom_distance = np.nan_to_num(Atom_distance)
            VdW_cluster[str(pdb_file).split('-')[1]] = cluster_calc(Atom_distance, by_atom)
            
            from helper_function import intersect_vol
            Atom_distance_nan = np.where(Atom_distance==0, np.nan, Atom_distance)
            Atom_volume = intersect_vol(Atom_distance_nan, 6, 6)
            VdW_volume[str(pdb_file).split('-')[1]] = Atom_volume
            
    return VdW_cluster, VdW_volume
                
                
def H_bond_calc(path, pqr_files=None):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    from helper_function import distance
    from helper_function import angle_calc
    
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

            
        from helper_function import distance
        from helper_function import angle_calc
        angle = angle_calc(Donor_array, H_array, Acceptor_array)
        HB_dict[str(pqr_file).split('.')[0]] = angle

    return HB_dict
                                            

def AA2s4pred (directory_S4pred, output_path, AA_seq, prot):
    import os
    from helper_function import fasparse
    os.getcwd()
    # call s4pred and create fas file
    fastapath = os.path.join(output_path,f'{prot}.fasta')
    faspath = os.path.join(output_path,f'{prot}.fas')
    abs_fasta = os.path.abspath(fastapath)
    abs_fas = os.path.abspath(faspath)

    if os.path.isfile(fastapath):
        print(f'fasta file already exists {fastapath}')
    else:
        with open(os.path.join(output_path,f'{prot}.fasta'), "w") as fasta_file:
            fasta_file.write(f">{prot}\n{AA_seq}\n")
        print(f'fasta file created {fastapath}')
        
    if os.path.isfile(faspath):
        print('fas file already exists')
    else:
        os.chdir(directory_S4pred)
        os.system(f'python3 run_model.py "{abs_fasta}" > "{abs_fas}"')
        os.chdir('../../')
        print('fas file created')
     
    #read fas file and   
    sec_pred = fasparse(abs_fas)    

    return sec_pred 











#-------------------------------------------------Mutation Creation-------------------------------------------------#             
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
    from function import salt_bridge
    from function import H_bond_calc
    from function import VdW_interaction
    from function import pdb2pqr
    from helper_function import remove_nan
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
    print('Salt_bridge finished')
    H_bond = H_bond_calc(output_path, pqr_file)
    print('H_bond finished')
    VdW_clust, VdW_vol = VdW_interaction(input_path, pdb_file, by_atom = True)
    print('VdW_interaction finished')
    
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
            
  




#function for predicting alpha helices, beta sheets and turns, turns not tested yet, helices and sheets testes on one protein, worked, more to come, improvements will follow
#depends on dicts:
helixvalues = {'E':1.59,'A':1.41,'L':1.34,'M':1.3,'Q':1.27,'K':1.23,'R':1.21,'H':1.05,'V':0.9,'I':1.09,'Y':0.74,'C':0.66,'W':1.02,'F':1.16,'T':0.76,'G':0.43,'N':0.76,'P':0.34,'S':0.57,'D':0.99,'U':0.66}
sheetvalues = {'E':0.52,'A':0.72,'L':1.22,'M':1.14,'Q':0.98,'K':0.69,'R':0.84,'H':0.8,'V':1.87,'I':1.67,'Y':1.45,'C':1.4,'W':1.35,'F':1.33,'T':1.17,'G':0.58,'N':0.48,'P':0.31,'S':0.96,'D':0.39,'U':1.4}
loopvalues = {'E':1.01,'A':0.82,'L':0.57,'M':0.52,'Q':0.84,'K':1.07,'R':0.9,'H':0.81,'V':0.41,'I':0.47,'Y':0.76,'C':0.54,'W':0.65,'F':0.59,'T':0.9,'G':1.77,'N':1.34,'P':1.32,'S':1.22,'D':1.24,'U':0.54}

def univt2(seq:str, size:float):
    helic=[]
    sheet=[]
    turn=[]
    for t in range(size,len(seq)-size):
        frame = [*seq[t - size:t + size+1]]                         #define neighbors of t
        mnh = sum(list(map(helixvalues.get,frame))) / (2*size+1)
        mns = sum(list(map(sheetvalues.get,frame))) / (2*size+1)
        mnt = sum(list(map(loopvalues.get,frame))) / (2*size+1)
        if mnh > 1.1 and mnh > mns and mnh > mnt:
            helic.append(t)
        elif mns > 1 and mns > mnh and mns > mnt:
            sheet.append(t)
        elif mnt > 1 and mnt > mnh and mnt > mns:
            turn.append(t)   
    counth=0
    counts=0
    countt=0
    for p in range(1,len(helic)-n):
        if helic[p+2]==helic[p]+2 and helic[p-1]!=helic[p]-1:
            counth+=1
    if len(helic) > 2:
        if helic[2] == helic[0]+2:
            counth+=1
    for n in range(1,len(sheet)-n):
        if sheet[n+2]==sheet[n]+2 and sheet[n-1]!=sheet[n]-1:
            counts+=1
    if len(sheet) >= 2:
        if sheet[2] == sheet[0]+2:
            counts+=1
    for s in range(1,len(turn)-n):
        if turn[s+2]==turn[s]+2 and turn[s-1]!=turn[s]-1:
            countt+=1
    if len(turn) > 2:
        if turn[2] == turn[0]+2:
            countt+=1
    
    return [counth,counts,countt]
    '''print(f'Helices:{counth}')
    print(f'Sheets:{counts}')
    print(f'Turns:{countt} (Nicht getestet)')'''
     
     
#fixed size of 2   
def univt3(seq:str):
    helic=[]
    sheet=[]
    turn=[]
    for t in range(2,len(seq)-2):
        frame = [*seq[t - 2:t + 3]]                         #define neighbors of t
        mnh = sum(list(map(helixvalues.get,frame))) / 5
        mns = sum(list(map(sheetvalues.get,frame))) / 5
        mnt = sum(list(map(loopvalues.get,frame))) / 5
        if mnh > 1.1 and mnh > mns and mnh > mnt:
            helic.append(t)
        elif mns > 1 and mns > mnh and mns > mnt:
            sheet.append(t)
        elif mnt > 1 and mnt > mnh and mnt > mns:
            turn.append(t) 
    counth=0
    counts=0
    countt=0
    for p in range(1,len(helic)-2):
        if helic[p+2]==helic[p]+2 and helic[p-1]!=helic[p]-1:
            counth+=1
    if len(helic) > 2:
        if helic[2] == helic[0]+2:
            counth+=1
    for n in range(1,len(sheet)-2):
        if sheet[n+2]==sheet[n]+2 and sheet[n-1]!=sheet[n]-1:
            counts+=1
    if len(sheet) >= 2:
        if sheet[2] == sheet[0]+2:
            counts+=1
    for s in range(1,len(turn)-2):
        if turn[s+2]==turn[s]+2 and turn[s-1]!=turn[s]-1:
            countt+=1
    if len(turn) > 2:
        if turn[2] == turn[0]+2:
            countt+=1
    
    return [counth,counts,countt]
    '''print(f'Helices:{counth}')
    print(f'Sheets:{counts}')
    print(f'Turns:{countt} (Nicht getestet)')'''
        
def p_val(corr, n, alpha):
    import math
    import scipy.stats as stats
    if math.sqrt((1-(corr**2))/(n-2)) != 0 and n-2 != 0:
        t = (corr)/(math.sqrt((1-(corr**2))/(n-2)))
        p = 1 - stats.t.cdf(t, n-2)
        return [p, p < alpha]

def pain(data, by= None, lower = 0.1, upper = 0.9):
    """
    data: pandas DataFrame or Series, if series then no argument by \n
    by: column to filter by (string) (only applicable if data is a DataFrame) \n
    upper and lower: threshold percentages, for example lower = 0.1, upper = 0.9"""
    import pandas as pd
    if isinstance(data, pd.DataFrame):
        lower_threshold = data[by].quantile(lower)
        upper_threshold = data[by].quantile(upper)
        return data[(data[by] <= lower_threshold) | (data[by] >= upper_threshold)].reset_index(drop=True)
    if isinstance(data, pd.Series):
        lower_threshold = data.quantile(lower)
        upper_threshold = data.quantile(upper)
        return data[(data <= lower_threshold) | (data >= upper_threshold)].reset_index(drop=True)
    else:
        raise ValueError('Maxi ist ein Idiot')
