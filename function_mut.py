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
    from function import rel_aa_comp
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
    from function_mut import mut_apply
    from function_mut import diff_weighted
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
        
        # Generate all combinations of substitutions
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
    from function import rel_aa_comp
    from function_mut import mut_apply
    from function_mut import mut_live_test
    from Aminoacid_lists import AA_polar_neutral
    
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
                    if key in helix:
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
                    if key in helix:
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