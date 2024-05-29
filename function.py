#function that calculates the relative amino acid composition based on sequence entry and definded amino acid property
def rel_aa_comp(Sequence:str, AA_property:str) -> list: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count/len(Sequence)

def rel_aa(Sequence:str, AA_property:str) -> list: 
    count = 0
    for element in Sequence:  
        if element in AA_property:
            count += 1
    return count

def salt_bridge(path, pdb_files=None):
    if pdb_files is None:
        pdb_files = [f for f in os.listdir(path) if f.endswith('.pdb')]
    if isinstance(pdb_files, str):
        pdb_files = [pdb_files] 
    for pdb_file in pdb_files:
        Salt_bridges = dict()
        Asp_Glu_array = np.empty((0, 4))
        Lys_Arg_His_array = np.empty((0, 4))
        
        with open(os.path.join(path, str(pdb_file))) as f:
            Salt= dict()
            Asp_Glu_array = np.empty((0,4))
            Lys_Arg_His_array = np.empty((0,4))
            for line in f:
                line = line.strip()
                if line.startswith('ATOM'):
                    if ('ASP' in line and 'OD' in line) or ('GLU' in line and 'OE' in line):
                        line_array = np.array([[line[8:11].strip(), line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]])
                        line_array = line_array.astype('float64')
                        Asp_Glu_array = np.append(Asp_Glu_array, line_array, axis = 0)
                    if ('LYS' in line and 'NZ' in line) or ('ARG' in line and 'NH' in line) or ('HIS' in line and 'NE' in line) or ('HIS' in line and 'ND' in line):
                        line_array = np.array([[line[8:11].strip(), line[27:38].strip(), line[39:46].strip(), line[47:54].strip()]])
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
