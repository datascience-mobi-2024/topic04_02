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

            #calculate distance clean up array
            distance = cdist(Asp_Glu_array[:,1:], Lys_Arg_His_array[:,1:], metric='euclidean') #calculate distance
            distance = np.concatenate((np.array([Lys_Arg_His_array[:,0]]), distance), axis=0) #add atom number from Lys_Arg_His to array
            distance = np.concatenate((np.insert(np.array([Asp_Glu_array[:,0]]), 0, None).reshape(-1,1), distance), axis=1) #add atom number from Asp_Glu to array
            distance[1:, 1:][distance[1:, 1:] > 4] = np.nan #set distance > 4 to nan
            rows_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=1)]),0, None) #find rows with all nan values
            cols_with_nan = np.insert(np.array([np.all(np.isnan(distance[1:, 1:]), axis=0)]),0, None) #find columns with all nan values
            distance = distance[~rows_with_nan, :] #delete rows with all nan values
            distance = distance[:, ~cols_with_nan] #delete columns with all nan values
            Salt_bridges[str(pdb_file).split('-')[1]] = distance

    return Salt_bridges
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
