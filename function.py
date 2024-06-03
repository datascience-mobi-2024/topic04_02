#function that calculates the relative amino acid composition based on sequence entry and definded amino acid property
def rel_aa_comp(Sequence:str, AA_property:str) -> str: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count/len(Sequence)

def rel_aa(Sequence:str, AA_property:str) -> str: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count

########################### Used for 3D structure analysis ###########################

#https://www.bioinformation.net/003/002800032008.pdf
def salt_bridge(path, pdb_files=None):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
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
                line = line.strip()
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

def VdW_interaction(path, pdb_files=None, output = None):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files]
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
            VdW_cluster = {}
            Atom_distance = np.nan_to_num(Atom_distance)
            VdW_cluster = cluster_calc(Atom_distance)
            
            from helper_function import intersect_vol
            VdW_volume = {}
            Atom_distance_nan = np.where(Atom_distance==0, np.nan, Atom_distance)
            Atom_volume = intersect_vol(Atom_distance_nan, 6, 6)
            VdW_volume[str(pdb_file).split('-')[1]] = Atom_volume
            
    return VdW_cluster, VdW_volume
                
                
def VdW_interaction(path, pqr_files=None, output = None):
    import numpy as np
    import os
    import scipy
    from scipy.spatial.distance import cdist
    Donor_list = [('GLN', 'NE2', 'HE21'), ('GLN', 'NE2', 'HE22'), ('GLU', 'OE2', 'HE2'), ('ASP', 'OD2', 'HD2'), ('ASN', 'ND2', 'HD21'), ('ASN', 'ND2', 'HD22'),
                    ('HIS', 'NE2', 'HE2'), ('HIS', 'ND1', 'HD2'), ('LYS', 'NZ', 'HZ1'), ('LYS', 'NZ', 'HZ2'), ('LYS', 'NZ', 'HZ3'), ('ARG', 'NE', 'HE'), 
                    ('ARG', 'NH1', 'HH11'), ('ARG', 'NH1', 'HH12'), ('ARG', 'NH2', 'HH21'), ('ARG', 'NH2', 'HH22'), ('SER', 'OG', 'HG'), ('THR', 'OG1', 'HG1'), 
                    ('TRP', 'NE1', 'HE1'), ('TYR', 'OH', 'HH')]
    Acceptor_list = [('GLN', 'OE1'), ('ASP', 'OD1'), ('ASP', 'OD2'), ('ASN', 'OD1')]
    
    if pqr_files is None:
        pqr_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pqr_files, str):
        pqr_files = [pqr_files]
    for pqr_file in pqr_files:
        with open(os.path.join(path, 'C0H3Z2.pqr')) as f:
            HB_dic = {}
            Donor_array = np.empty((0, 4))
            H_array = np.empty((0, 4))
            Acceptor_array = np.empty((0, 4))
            aa_cache = []
            atom_cache = []
            test_count = 0
            for line in f:
                if line.startswith('ATOM'):
                    if not aa_cache:
                        aa_cache.append(line.split()[3])
                        aa_cache.append(line.split()[4])
                        atom_cache.append(line)
                    elif aa_cache[1] == line.split()[4]:
                        atom_cache.append(line)
                    elif aa_cache[1] != line.split()[4]:
                        test_count += 1             
                        for n in range(len(Donor_list)):
                            if aa_cache[0] == Donor_list[n][0]:
                                for i in range(len(atom_cache)):
                                    if Donor_list[n][1] == atom_cache[i].split()[2]:
                                        if ('GLU' == aa_cache[0]) and 'OE2' in atom_cache[i].split()[2]:
                                            if any('HE2' in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Donor_array = np.append(Donor_array, line_array, axis=0) 
                                        elif ('ASP' == aa_cache[0]) and 'OD2' in atom_cache[i].split()[2]:
                                            if any('HD2' in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Donor_array = np.append(Donor_array, line_array, axis=0)                                             
                                        elif ('HIS' == aa_cache[0]) and 'ND1' in atom_cache[i].split()[2]:
                                            if any ('HD1' in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Donor_array = np.append(Donor_array, line_array, axis=0)
                                        elif ('HIS' == aa_cache[0]) and 'NE2' in atom_cache[i].split()[2]:
                                            if any('HE2' in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Donor_array = np.append(Donor_array, line_array, axis=0)
                                        else: 
                                            line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                            line_array = line_array.astype('float64')
                                            Donor_array = np.append(Donor_array, line_array, axis=0)         
                                    elif Donor_list[n][2] == atom_cache[i].split()[2]:
                                        line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                        line_array = line_array.astype('float64')
                                        H_array = np.append(H_array, line_array, axis=0)
                        for n in range(len(Acceptor_list)):
                            if aa_cache[0] == Acceptor_list[n][0]:
                                for i in range(len(atom_cache)):
                                    if Acceptor_list[n][1] == atom_cache[i].split()[2]:
                                        if ('GLU' == aa_cache[0]) and 'OE2' in atom_cache[i].split()[2]:
                                            if any('HE2' not in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Acceptor_array = np.append(Acceptor_array, line_array, axis=0)
                                        elif ('ASP' == aa_cache[0]) and 'OD' in atom_cache[i].split()[2]:
                                            if any('HD' not in string for string in atom_cache):
                                                line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                                line_array = line_array.astype('float64')
                                                Acceptor_array = np.append(Acceptor_array, line_array, axis=0)
                                        else: 
                                            line_array = np.array([[atom_cache[i].split()[1], atom_cache[i].split()[5], atom_cache[i].split()[6], atom_cache[i].split()[7]]])
                                            line_array = line_array.astype('float64')
                                            Acceptor_array = np.append(Acceptor_array, line_array, axis=0) 
            from helper_function import distance
            from helper_function import angle_calc
            angle = angle_calc(Donor_array, H_array, Acceptor_array)
            HB_dic[str(pqr_file).split('.')[1]] = angle
    return HB_dic
                                            

             
                
            
  




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
